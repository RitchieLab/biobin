/*
 * PopulationManager.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#include "PopulationManager.h"

#include "binapplication.h"
#include "binmanager.h"

#include "tests/Test.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <math.h>

#include <iostream>
#include <limits>

using std::fill;
using std::vector;
using std::string;
using std::map;
using std::ostream;
using std::ifstream;
//using boost::array;
using boost::unordered_map;
using boost::dynamic_bitset;

using Knowledge::Locus;
using Knowledge::Region;
using Knowledge::Information;

using BioBin::Test::Test;
using BioBin::Utility::Phenotype;

using boost::program_options::validation_error;
using boost::algorithm::split;
using boost::algorithm::is_any_of;
using boost::trim;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using Knowledge::Locus;

namespace BioBin{

float PopulationManager::c_phenotype_control = 0;
string PopulationManager::c_phenotype_file = "";
string PopulationManager::c_covariate_file = "";
float PopulationManager::c_min_control_frac = 0.125;
PopulationManager::DiseaseModel PopulationManager::c_model =
		PopulationManager::ADDITIVE;
PopulationManager::WeightModel PopulationManager::c_weight_type =
		PopulationManager::MAX;
vector<BioBin::Test::Test*> PopulationManager::c_tests;

bool PopulationManager::RareCaseControl = true;
bool PopulationManager::c_use_calc_weight = false;
bool PopulationManager::NoSummary = false;

PopulationManager::PopulationManager(const string& vcf_fn) :
		_vcf_fn(vcf_fn), _use_custom_weight(false), _info(0){
}

unsigned int PopulationManager::genotypeContribution(const Locus& loc) const{
	unordered_map<const Locus*, bitset_pair>::const_iterator itr = _genotypes.find(&loc);

	if(itr != _genotypes.end()){
		return getTotalContrib((*itr).second);
	} else {
		return 0;
	}
}

float PopulationManager::getAvgGenotype(const Locus& loc, const dynamic_bitset<>* nonmiss_status) const{
	unordered_map<const Locus*, bitset_pair>::const_iterator itr = _genotypes.find(&loc);
	unsigned int n_vars = 0;
	float nonmiss = 1;

	if(itr != _genotypes.end()){
		n_vars =  getTotalContrib((*itr).second, nonmiss_status);
		if(nonmiss_status){
			nonmiss = (*nonmiss_status & (~((*itr).second.first & (*itr).second.second))).count();
		} else {
			nonmiss = (~((*itr).second.first & (*itr).second.second)).count();
		}
	}

	// will return 0 for a missing locus
	return n_vars / nonmiss;
}

unsigned short PopulationManager::getIndivGeno(const Locus& loc, int pos) const{
	unordered_map<const Knowledge::Locus*, bitset_pair>::const_iterator it = _genotypes.find(&loc);
	if(it == _genotypes.end()){
		return missing_geno;
	}

	unsigned short n_var = 2*(*it).second.first[pos] + (*it).second.second[pos];

	if(n_var == 3){
		// if we're here, then we had both bits set to 1, which means missing
		n_var = missing_geno;
	}else if(c_model == DOMINANT){
		n_var = n_var > 0;
	} else if (c_model == RECESSIVE){
		n_var = n_var > 1;
	}

	return n_var;

}

bool PopulationManager::isRare(const Locus& locus, const bitset_pair& status, float lower, float upper) const{
	bool rare = false;
	float currmaf;

	unordered_map<const Locus*, bitset_pair>::const_iterator g_itr = _genotypes.find(&locus);
	if(g_itr != _genotypes.end()){
		currmaf = getMAF((*g_itr).second, & status.first);
		rare = currmaf <= upper && currmaf >= lower;

		if(!rare && RareCaseControl){
			currmaf = getMAF((*g_itr).second, & status.second);
			rare = currmaf <= upper && currmaf >= lower;
		}
	}

	return rare;
}

float PopulationManager::getPhenotypeVal(const string& sample, const Phenotype& pheno) const{
	static const float missing_status = std::numeric_limits<float>::quiet_NaN();
	boost::unordered_map<std::string, std::vector<float> >::const_iterator pheno_status = _phenos.find(sample);
	float status = missing_status;
	if (pheno_status != _phenos.end()){
		status = (*pheno_status).second[pheno.getIndex()];
	}
	return status;

}

unsigned int PopulationManager::getSamplePosition(const std::string& s) const{
	unsigned int pos = static_cast<unsigned int>(-1);
	boost::unordered_map<std::string, unsigned int>::const_iterator pos_itr = _positions.find(s);
	if (pos_itr != _positions.end()){
		pos = (*pos_itr).second;
	}
	return pos;
}

const vector<float>& PopulationManager::getCovariates(const string& sample) const{
	static const vector<float> missing_covars;
	boost::unordered_map<std::string, std::vector<float> >::const_iterator covar_status = _covars.find(sample);

	return ( (covar_status == _covars.end()) ? missing_covars : (*covar_status).second );

}

float PopulationManager::getTotalIndivContrib(const Bin& b,int pos, const Phenotype& pheno) const {
	// Accumulate the contribution of this person in this bin
	Bin::const_locus_iterator l_itr = b.variantBegin();
	float bin_count = 0;
	while(l_itr != b.variantEnd()){
		bin_count += getIndivContrib(**l_itr, pos, pheno, true, b.getRegion());
		++l_itr;
	}
	return bin_count;
}

unsigned int PopulationManager::readVCFHeader(std::istream& v){
	unsigned int lineno = 0;
	bool header_read = false;
	string curr_line;

	while (!header_read && getline(v, curr_line)) {
		++lineno;
		vector<string> fields;
		if (curr_line.size() > 2 && curr_line[0] == '#' && curr_line[1] != '#') {
			// In this case, we are looking at the "#CHROM" line (we hope!)
			split(fields, curr_line, boost::is_any_of("\t"));
			unsigned int n_fields = fields.size();

			// Sanity check here! Note that I never use the "QUAL" or "INFO"
			// fields, so I'm not going to check for them!
			if (n_fields < 9 || fields[0] != "#CHROM" || fields[1] != "POS"
					|| fields[2] != "ID" || fields[3] != "REF" || fields[4]
					!= "ALT" || fields[6] != "FILTER" || fields[8] != "FORMAT") {
				std::cerr << "ERROR: VCF header line is malformed, please check your VCF input."
						<< std::endl;

				//throw exception here
				throw std::runtime_error("VCF Header malformed");
			}

			_sample_names.insert(_sample_names.end(), fields.begin() + 9, fields.end());
			for(unsigned int i=0; i<fields.size()-9; i++){
				_positions.insert(std::make_pair(fields[i+9], i));
			}

			header_read = true;
		}
	}

	return lineno;
}

void PopulationManager::loadIndividuals(){

	if(c_covariate_file != ""){
		parseTraitFile(c_covariate_file, _covar_names, _covars, "covar");
	}

	bool allControl = false;
	if(c_phenotype_file != ""){
		parseTraitFile(c_phenotype_file, _pheno_names, _phenos, "pheno");
	} else {
		allControl = true;
	}

	// Now, we go through and determine who is a case and who is a control
	boost::unordered_map<string, unsigned int>::const_iterator p_itr = _positions.begin();
	boost::unordered_map<string, unsigned int>::const_iterator p_end = _positions.end();
	boost::unordered_map<string, vector<float> >::const_iterator pheno_itr;
	boost::unordered_map<string, vector<float> >::const_iterator pheno_not_found = _phenos.end();
	if(allControl){
		_pheno_names.push_back("");
		// this will create the vector and make it length 1 (w/ value 0) all at once!
		while(p_itr != p_end){
			_phenos[(*p_itr).first].push_back(0);
			++p_itr;
		}
		dynamic_bitset<> cab = dynamic_bitset<>(_positions.size());
		dynamic_bitset<> cob = dynamic_bitset<>(_positions.size());
		cob.set();
		_pheno_status.push_back(std::make_pair(cob, cab));
	} else {

		// push back all phenotypes that consists of "all missing"
		dynamic_bitset<> cab(_positions.size());
		dynamic_bitset<> cob(_positions.size());
		_pheno_status.reserve(_pheno_names.size());
		for(unsigned int i=0; i<_pheno_names.size(); i++){
			_pheno_status.push_back(std::make_pair(cob, cab));
		}

		// now, iterate over ALL people that we have phenotypes for
		while(p_itr != p_end){
			pheno_itr = _phenos.find((*p_itr).first);

			if (pheno_itr != pheno_not_found){
				for(unsigned int i=0; i<std::min((*pheno_itr).second.size(), _pheno_names.size()); i++){
					if((*pheno_itr).second[i] == c_phenotype_control){
						_pheno_status[i].first.set((*p_itr).second);
					}else{
						_pheno_status[i].second[(*p_itr).second] = !std::isnan((*pheno_itr).second[i]);
					}
				}
			}
			++p_itr;
		}
	}

	// now, check each phenotype to make sure we have enough controls!
	for(unsigned int i=0; i<_pheno_names.size(); i++){
		unsigned int n_controls = _pheno_status[i].first.count();
		unsigned int n_cases = _pheno_status[i].second.count();
		unsigned int total = n_controls + n_cases;

		// If we don't have enough controls, print a warning and make everyone a control
		if ((n_controls / (float) total) < c_min_control_frac){
			std::cerr << "WARNING: In Phenotype '" << _pheno_names[i]
			          << "', number of controls is less than "
			          << c_min_control_frac * 100
			          << "% of the data.  Using all individuals as controls"
			          << std::endl;
			_pheno_status[i].first |= _pheno_status[i].second;
			_pheno_status[i].second.reset();
			n_controls = total;
			n_cases = 0;
		}else if(1-(n_controls / (float) total) < c_min_control_frac && (n_controls / total != 1)){
			std::cerr << "WARNING: In phenotype '" << _pheno_names[i]
			          << "', number of cases is less than "
			          << c_min_control_frac * 100
			          << "% of the data.  Allele frequencies for cases may be unreliable"
			          << std::endl;
		}

		// Print a warning if rare variants will only be fixed variants
		if (1 / static_cast<float>(2 *n_controls) > BinManager::mafCutoff){
			std::cerr << "WARNING: MAF cutoff is set so low that only variants "
					  << "fixed in controls are rare for phenotype '"
					  << _pheno_names[i] << "'." << std::endl;
		}
		if (RareCaseControl && n_controls != total && 1 / static_cast<float>(2*(n_cases)) > BinManager::mafCutoff){
			std::cerr << "WARNING: MAF cutoff is set so low that only variants "
					  << "fixed in cases are rare for phenotype '"
					  << _pheno_names[i] << "." << std::endl;
		}
	}
}

void PopulationManager::parseTraitFile(const string& filename,
		vector<string>& names_out,
		unordered_map<string, vector<float> >& vals_out,
		const string& var_prefix){

	// Open the file
	ifstream data_file(filename.c_str());
	if (!data_file.is_open()){
		std::cerr<<"WARNING: cannot find " << filename <<", ignoring.";
		return;
	}

	// Read the definition of the meta-group
	string line;
	vector<string> result;
	int line_no = 0;
	bool header_read = false;
	while(data_file.good()){
		getline(data_file, line);
		++line_no;
		trim(line);

		if(line.size() > 0){

			if(!header_read && line[0] == '#'){
				split(result, line, is_any_of(" \n\t"), boost::token_compress_on);
				for(unsigned int i=1; i<result.size(); i++){
					names_out.push_back(result[i]);
				}
			} else if(line[0] != '#'){
				split(result, line, is_any_of(" \n\t"), boost::token_compress_on);

				// if we haven't read the header, just assign arbitrary phenotype names
				if(!header_read){
					if(result.size() > 2){
						std::cerr << "WARNING: No header given for multiple traits, "
								<< "assigning sequential names" << std::endl;
						for(unsigned int i=1; i<result.size(); i++){
							names_out.push_back(var_prefix + "_" + boost::lexical_cast<string>(i));
						}
					}else{
						names_out.push_back("");
					}
				}

				if (result.size() != names_out.size() + 1){
					std::cerr << "WARNING: improperly formatted trait file " << filename
							  << " on line " << line_no << ", ignoring." << std::endl;
				}else if (_positions.find(result[0]) == _positions.end()){
						std::cerr << "WARNING: cannot find " << result[0] << " in VCF file." << std::endl;
				} else {
					vals_out[result[0]].reserve(names_out.size());
					for(unsigned int i=1; i<result.size(); i++){
						try{
							vals_out[result[0]].push_back(lexical_cast<float>(result[i]));
						}catch(bad_lexical_cast&){
							vals_out[result[0]].push_back(std::numeric_limits<float>::quiet_NaN());
						}
					}
				}
			}
			header_read = true;
		}
	}
}

float PopulationManager::getIndivContrib(const Locus& loc, int pos, const Utility::Phenotype& pheno, bool useWeights, const Region* const reg) const{

	unordered_map<const Knowledge::Locus*, bitset_pair>::const_iterator it = _genotypes.find(&loc);

	unsigned short n_var = getIndivGeno(loc, pos);

	float wt = 1;
	if(n_var == missing_geno){
		n_var = 0;
	}

	// Cache the weights so we aren't wasting so much effort.
	// Also, only calculate the weight if n_var > 0 (o/w will multiply out to 0)
	if(n_var != 0 && useWeights){
		wt = getLocusWeight(loc, pheno, reg);
	}

	return n_var * wt;
}

unsigned int PopulationManager::getTotalContrib(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss)  const{

	boost::dynamic_bitset<> nonmissing = ~(geno.first & geno.second);
	if(nonmiss != 0){
		nonmissing &= *nonmiss;
	}

	switch (c_model) {
	case ADDITIVE:
		return 2 * (nonmissing & geno.first).count() + (nonmissing
				& geno.second).count();
	case DOMINANT:
		return (nonmissing & geno.first).count() + (nonmissing
				& geno.second).count();
	case RECESSIVE:
		return (nonmissing & geno.first).count();
	default:
		return 0;
	}
}

float PopulationManager::getMAF(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss) const{
	boost::dynamic_bitset<> nonmissing = ~(geno.first & geno.second);
	if(nonmiss != 0){
		nonmissing &= *nonmiss;
	}

	float maf = (2 * (nonmissing & geno.first).count() + (nonmissing & geno.second).count()) / static_cast<float>(2*nonmissing.count());

	return std::min(maf, 1-maf);
}

float PopulationManager::calcBrowningWeight(unsigned long N, unsigned long M) const{
	// Madsen + Browning Weight
	// A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic
	// PLOSGenetics, Feb 2009, Vol. 5, Issue 2, e10000384

	return (2*N+2)/sqrt(N*(M+1)*(2*N-M+1));

	// If given the frequency as opposed to the number of mutations:
	//return 2*sqrt((1+N)*(1+N)/(N+2*N*N+4*F*(1-F)*N*N*N));
}

float PopulationManager::getLocusWeight(const Locus& loc, const Phenotype& pheno, const Region* reg) const{

	float wt = 1;

	if(_use_custom_weight){
		wt *= getCustomWeight(loc, reg);
	}

	if(c_use_calc_weight){
		wt *= calcWeight(loc, pheno.getStatus());
	}

	return wt;

}

float PopulationManager::calcWeight(const Locus& loc, const bitset_pair& status) const{
	// *_a = Affected (cases), *_u = Unaffected (controls)
	// N_* = Number (population), F_* = Frequency
	const boost::dynamic_bitset<>& control_bitset=status.first;
	const boost::dynamic_bitset<>& case_bitset=status.second;

	int N_a, N_u, N;
	int M_a, M_u, M;
	float w_a = std::numeric_limits<float>::quiet_NaN();
	float w_u = 1;

	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotypes.find(&loc);

	if(it == _genotypes.end()){
		return 1;
	}

	dynamic_bitset<> nonmissing = ~((*it).second.first & (*it).second.second);

	N_u = (nonmissing & control_bitset).count() ;
	N_a = (nonmissing & case_bitset).count();
	N = N_u + N_a;

	M_u = 2*((*it).second.first & control_bitset & nonmissing).count() +
			((*it).second.second & control_bitset & nonmissing).count();
	M_a = 2*((*it).second.first & case_bitset & nonmissing).count() +
			((*it).second.second & case_bitset & nonmissing).count();
	M = M_u + M_a;

	float weight = 1;
	if(c_weight_type != OVERALL){
		if(N_u){
			w_u = calcBrowningWeight(N_u, M_u);
		}
		if(c_weight_type != CONTROL && N_a){
			w_a = calcBrowningWeight(N_a, M_a);
		}

		switch(c_weight_type){
		// NOTE: since w_a may be NaN, and all comparisons with NaN are by
		// definition false, thenw_a MUST be returned ONLY when the expression
		// evaluates to true, because that means that w_a CANNOT be NaN in that case
		case MAX:
			weight = (N_u > 0 && w_u < w_a) || (N_u == 0 && w_a == w_a) ? w_a : w_u;
			break;
		case MIN:
			weight = (N_u > 0 && w_u > w_a) || (N_u == 0 && w_a == w_a) ? w_a : w_u;
			break;
		case CONTROL:
			weight = w_u;
			break;
		default:
			// WE SHOULD NOT BE HERE
			weight = 1;
		}
	} else if(N){ // c_weight_type == OVERALL
		weight = calcBrowningWeight(N, M);
	}

	return weight;
}

float PopulationManager::getCustomWeight(const Locus& loc, const Region* const reg) const{
	static float weight_cache = 1;
	static const Locus* loc_cache = 0;
	static const Region* reg_cache = 0;


	if(_info){
		if(&loc != loc_cache || reg != reg_cache){
			loc_cache = &loc;
			reg_cache = reg;
			weight_cache = _info->getSNPWeight(loc, reg);
		}
	}
	return weight_cache;
}

boost::array<unsigned int, 2> PopulationManager::getBinCapacity(Bin& bin, const bitset_pair& status) const {

	const boost::dynamic_bitset<>& control_bitset=status.first;
	const boost::dynamic_bitset<>& case_bitset=status.second;
	Bin::const_locus_iterator b_itr = bin.variantBegin();
	Bin::const_locus_iterator b_end = bin.variantEnd();
	boost::array<unsigned int, 2> capacity;
	capacity[0] = 0;
	capacity[1] = 0;

	unordered_map<const Locus*, bitset_pair>::const_iterator l_end = _genotypes.end();
	unordered_map<const Locus*, bitset_pair>::const_iterator l_itr;
	while (b_itr != b_end) {

		l_itr = _genotypes.find(*b_itr);
		if (l_itr != l_end) {
			dynamic_bitset<> nonmiss = ~((*l_itr).second.first & (*l_itr).second.second);
			capacity[0] += (nonmiss & control_bitset).count();
			capacity[1] += (nonmiss & case_bitset).count();
		}
		++b_itr;
	}

	capacity[0] *= (1 + (c_model == ADDITIVE));
	capacity[1] *= (1 + (c_model == ADDITIVE));

	return capacity;
}

void PopulationManager::printEscapedString(ostream& os, const string& toPrint, const string& toRepl, const string& replStr) const{
	os << boost::algorithm::replace_all_copy(toPrint, toRepl, replStr);
}

string PopulationManager::getEscapeString(const string& sep) const{
	string sep_repl = "_";
	if (sep == sep_repl){
		sep_repl = "-";
	}
	return sep_repl;
}

void PopulationManager::printBinsTranspose(std::ostream& os, const BinManager& bins, const Phenotype& pheno, const std::string& sep) const{
	static const float missing_status = std::numeric_limits<float>::quiet_NaN();

	std::string sep_repl = getEscapeString(sep);

	// Print 1st line
	printEscapedString(os, "Bin Name", sep, sep_repl);
	os << sep;
	if (!NoSummary) {
		printEscapedString(os, "Total Variants", sep, sep_repl);
		os << sep;
		printEscapedString(os, "Total Loci", sep, sep_repl);
		os << sep;
		printEscapedString(os, "Control Loci Totals", sep, sep_repl);
		os << sep;
		printEscapedString(os, "Case Loci Totals", sep, sep_repl);
		os << sep;
		printEscapedString(os, "Control Bin Capacity", sep, sep_repl);
		os << sep;
		printEscapedString(os, "Case Bin Capacity", sep, sep_repl);
	}

	for(unsigned int i=0; i<c_tests.size(); i++){
		os << sep;
		printEscapedString(os, c_tests[i]->getName() + " p-value", sep, sep_repl);
	}

	boost::unordered_map<std::string, unsigned int>::const_iterator m_itr = _positions.begin();
	while(m_itr != _positions.end()){
		os << sep;
		printEscapedString(os, (*m_itr).first, sep, sep_repl);
		++m_itr;
	}

	os << "\n";

	printEscapedString(os, "Status", sep, sep_repl);
	os << sep;
	if(!NoSummary){
		os << missing_status << sep << missing_status << sep << missing_status
		   << sep << missing_status << sep << missing_status << sep << missing_status;
	}

	m_itr = _positions.begin();
	while(m_itr != _positions.end()){
		os << sep << getPhenotypeVal((*m_itr).first, pheno);
		++m_itr;
	}

	os << "\n";

	BinManager::const_iterator b_itr = bins.begin();
	Bin::const_locus_iterator l_itr;
	boost::unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator loc_itr;

	vector<vector<double> > test_pvals(c_tests.size());
	for(unsigned int i=0; i<c_tests.size(); i++){
		test_pvals[i].reserve(bins.size());
		c_tests[i]->runAllTests(*this, pheno,bins,test_pvals[i]);
	}

	unsigned int i=0;
	while (b_itr != bins.end()) {
		// Print bin name
		printEscapedString(os, (*b_itr)->getName(), sep, sep_repl);

		if (!NoSummary) {
			// print total var
			os << sep << (*b_itr)->getSize();

			// print total loci
			os << sep << (*b_itr)->getVariantSize();

			// print case/control loci
			boost::array<unsigned int, 2> num_loci;
			num_loci[0] = 0;
			num_loci[1] = 0;
			l_itr = (*b_itr)->variantBegin();
			while (l_itr != (*b_itr)->variantEnd()) {
				loc_itr = _genotypes.find((*l_itr));
				if (loc_itr != _genotypes.end()) {
					num_loci[0] += getTotalContrib((*loc_itr).second, & pheno.getStatus().first);
					num_loci[1] += getTotalContrib((*loc_itr).second, & pheno.getStatus().second);
				}
				++l_itr;
			}

			os << sep << num_loci[0] << sep << num_loci[1];

			// print case/control capacity
			boost::array<unsigned int, 2> capacity = getBinCapacity(**b_itr, pheno.getStatus());
			os << sep << capacity[0] << sep << capacity[1];
		} // End summary info

		// start test info
		for(unsigned int j=0; j<c_tests.size(); j++){
			os << sep;
			os << test_pvals[j][i];
		}


		// print for each person
		m_itr = _positions.begin();
		while (m_itr != _positions.end()) {
			os << sep << std::setprecision(4) << getTotalIndivContrib(**b_itr, (*m_itr).second, pheno);
			++m_itr;
		}
		os << "\n";
		++b_itr;
		++i;
	}
}

void PopulationManager::printBins(std::ostream& os, const BinManager& bins, const Phenotype& pheno, const std::string& sep) const{
	std::string sep_repl = getEscapeString(sep);
	static const float missing_status = std::numeric_limits<float>::quiet_NaN();

	BinManager::const_iterator b_itr = bins.begin();
	BinManager::const_iterator b_end = bins.end();

	// Print first line
	printEscapedString(os, "ID", sep, sep_repl);
	os << sep;
	printEscapedString(os, "Status", sep, sep_repl);
	while(b_itr != b_end){
		os << sep;
		printEscapedString(os, (*b_itr)->getName(), sep, sep_repl);
		++b_itr;
	}
	os << "\n";

	Bin::const_locus_iterator l_itr;
	Bin::const_locus_iterator l_end;

	boost::unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator loc_itr;
	boost::unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator loc_not_found =
			_genotypes.end();

	if(!NoSummary){

		// Print second Line (totals)
		printEscapedString(os, "Total Variants", sep, sep_repl);
		os << sep << missing_status;
		b_itr = bins.begin();
		b_end = bins.end();
		while(b_itr != b_end){
			os << sep << (*b_itr)->getSize();
			++b_itr;
		}
		os << "\n";

		// Print third line (variant totals)
		printEscapedString(os, "Total Loci", sep, sep_repl);
		os << sep << missing_status;
		b_itr = bins.begin();
		b_end = bins.end();
		while(b_itr != b_end){
			os << sep << (*b_itr)->getVariantSize();
			++b_itr;
		}
		os << "\n";

		int locus_count = 0;
		// Print 4th + 5th lines (variant totals for cases + controls
		for(int i=0; i<2; i++){
			printEscapedString(os, std::string(i ? "Case" : "Control") + " Loci Totals", sep, sep_repl);
			os << sep << missing_status;
			b_itr = bins.begin();
			b_end = bins.end();
			while(b_itr != b_end){
				l_itr = (*b_itr)->variantBegin();
				l_end = (*b_itr)->variantEnd();
				locus_count = 0;
				while(l_itr != l_end){
					loc_itr = _genotypes.find((*l_itr));
					if (loc_itr != loc_not_found){
						locus_count += getTotalContrib((*loc_itr).second, &(i ? pheno.getStatus().second : pheno.getStatus().first));
					}
					++l_itr;
				}
				os << sep << locus_count;
				++b_itr;
			}
			os << "\n";
		}

		// Print 6th + 7th Lines (bin capacities for cases and controls)
		for(int i=0; i<2; i++){
			printEscapedString(os, std::string(i ? "Case" : "Control") + " Bin Capacity", sep, sep_repl);
			os << sep << missing_status;
			b_itr = bins.begin();
			b_end = bins.end();
			while(b_itr != b_end){
				os << sep << getBinCapacity(**b_itr, pheno.getStatus())[i];
				++b_itr;
			}
			os << "\n";
		}
	}

	// Print the results of the tests
	for(unsigned int i=0; i<c_tests.size(); i++){
		printEscapedString(os, c_tests[i]->getName() + " p-value", sep, sep_repl);
		os << sep << missing_status;
		vector<double> pvals(bins.size());
		Test::Test* t = c_tests[i]->clone();
		t->runAllTests(*this, pheno,bins,pvals);
		delete t;
		for(unsigned int j=0; j<pvals.size(); j++){
			os << sep << pvals[j];
		}
		os << "\n";
	}

	boost::unordered_map<std::string, unsigned int>::const_iterator m_itr = _positions.begin();
	boost::unordered_map<std::string, unsigned int>::const_iterator m_end = _positions.end();

	int pos;
	float status;

	while (m_itr != m_end){
		b_itr = bins.begin();
		b_end = bins.end();

		pos = (*m_itr).second;
		status = getPhenotypeVal((*m_itr).first, pheno);

		printEscapedString(os, (*m_itr).first, sep, sep_repl);

		os << sep << status;

		while(b_itr != b_end){
			os << sep << std::setprecision(4) << getTotalIndivContrib(**b_itr, pos, pheno);
			++b_itr;
		}

		os << "\n";
		++m_itr;
	}
}





}

namespace std{
std::istream& operator>>(std::istream& in, BioBin::PopulationManager::DiseaseModel& model_out)
{
    std::string token;
    in >> token;
    if(token.size() > 0){
    	char s = token[0];
    	if(s == 'a' || s == 'A'){
    		model_out = BioBin::PopulationManager::ADDITIVE;
    	}else if(s == 'd' || s == 'D'){
    		model_out = BioBin::PopulationManager::DOMINANT;
    	}else if(s == 'r' || s == 'R'){
    		model_out = BioBin::PopulationManager::RECESSIVE;
    	}else{
    		throw validation_error(validation_error::invalid_option_value);
    	}
    }else{
    	throw validation_error(validation_error::invalid_option_value);
    }
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

ostream& operator<<(ostream& o, const BioBin::PopulationManager::DiseaseModel& m){
	o << (const char*) m;
	return o;
}

std::istream& operator>>(std::istream& in, BioBin::PopulationManager::WeightModel& model_out)
{
    std::string token;
    in >> token;
    if(token.size() > 0){
    	char s = token[0];
    	if(s == 'm' || s == 'M'){
    		if(token.size() > 1){
    			char s1 = token[1];
    			if(s1 == 'i' || s1 == 'I'){
    				model_out = BioBin::PopulationManager::MIN;
    			}else if (s1 == 'a' || s1 == 'A'){
    				model_out = BioBin::PopulationManager::MAX;
    			}else{
    				throw validation_error(validation_error::invalid_option_value);
    			}
    		}else{
    			throw validation_error(validation_error::invalid_option_value);
    		}
    	}else if(s == 'c' || s == 'C'){
    		model_out = BioBin::PopulationManager::CONTROL;
    	}else if(s == 'o' || s == 'O'){
    		model_out = BioBin::PopulationManager::OVERALL;
    	}else{
    		throw validation_error(validation_error::invalid_option_value);
    	}
    }else{
    	throw validation_error(validation_error::invalid_option_value);
    }
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

ostream& operator<<(ostream& o, const BioBin::PopulationManager::WeightModel& m){
	o << (const char*) m;
	return o;
}
}

