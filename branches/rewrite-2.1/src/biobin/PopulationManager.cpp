/*
 * PopulationManager.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#include "PopulationManager.h"

#include "binapplication.h"

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

using boost::program_options::validation_error;
using boost::algorithm::split;
using boost::algorithm::is_any_of;
using boost::trim;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using Knowledge::Locus;

namespace BioBin{

float PopulationManager::c_phenotype_control = 0;
vector<string> PopulationManager::c_phenotype_files;
float PopulationManager::c_min_control_frac = 0.125;
PopulationManager::DiseaseModel PopulationManager::c_model =
		PopulationManager::ADDITIVE;
PopulationManager::WeightModel PopulationManager::c_weight_type =
		PopulationManager::MAX;

//bool PopulationManager::CompressedVCF = false;
//bool PopulationManager::KeepCommonLoci = true;
bool PopulationManager::RareCaseControl = true;
//bool PopulationManager::OverallMajorAllele = true;
bool PopulationManager::c_use_calc_weight = false;
bool PopulationManager::_use_custom_weight = false;
bool PopulationManager::NoSummary = false;

PopulationManager::PopulationManager(const string& vcf_fn) : _vcf_fn(vcf_fn){
	//if(BinApplication::s_run_normal){
	//	vcf = new VCF::vcf_file(vcf_fn, CompressedVCF);
	//}
}

bool PopulationManager::isRare(const Locus& locus, float lower, float upper) const{
	bool rare = false;
	float currmaf;

	unordered_map<const Locus*, bitset_pair>::const_iterator g_itr = _genotypes.find(&locus);
	if(g_itr != _genotypes.end()){
		currmaf = getMAF((*g_itr).second, &_control_bitset);
		rare = currmaf <= upper && currmaf >= lower;

		if(!rare && RareCaseControl){
			currmaf = getMAF((*g_itr).second, &_case_bitset);
			rare = currmaf <= upper && currmaf >= lower;
		}
	}

	return rare;
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

			for(unsigned int i=0; i<fields.size()-9; i++){
				_positions.insert(std::make_pair(fields[i+9], i));
			}

			header_read = true;
		}
	}

	return lineno;
}

void PopulationManager::loadIndividuals(){

	// By default, everyone is missing who is not found in a phenotype file
	//_is_control = vector<bool>(size, true);
	_control_bitset = dynamic_bitset<>(_positions.size());
	_case_bitset = dynamic_bitset<>(_positions.size());

	// OK, now iterate through the phenotype files and load them up
	vector<string>::const_iterator itr = c_phenotype_files.begin();
	vector<string>::const_iterator end = c_phenotype_files.end();

	// if we give NO phenotype file, just set everyone as a control
	bool allControl = (itr == end);

	while(itr != end){
		parsePhenotypeFile(*itr);
		++itr;
	}

	// Now, we go through and determine who is a case and who is a control
	boost::unordered_map<string, int>::const_iterator p_itr = _positions.begin();
	boost::unordered_map<string, int>::const_iterator p_end = _positions.end();
	boost::unordered_map<string, float>::const_iterator pheno_itr;
	boost::unordered_map<string, float>::const_iterator pheno_not_found = _phenotypes.end();
	n_controls = 0;
	n_cases = 0;
	while(p_itr != p_end){
		pheno_itr = _phenotypes.find((*p_itr).first);

		if (pheno_itr != pheno_not_found){
			if((*pheno_itr).second == c_phenotype_control){

				++n_controls;
				_control_bitset.set((*p_itr).second);
			}else{
				++n_cases;
				_case_bitset[(*p_itr).second] = !std::isnan((*pheno_itr).second);
			}
		}
		++p_itr;
	}

	if(allControl){
		_control_bitset.set();
		_case_bitset.reset();
		n_controls = _control_bitset.size();
		n_cases = 0;
	}

	unsigned int total = n_controls + n_cases;

	// If we don't have enough controls, print a warning and make everyone a control
	if ((n_controls / (float) total) < c_min_control_frac){
		std::cerr << "WARNING: Number of controls is less than " <<
				c_min_control_frac * 100 << "% of the data.  Using all individuals as controls\n";
		_control_bitset |= _case_bitset;
		_case_bitset.reset();
	}else if(1-(n_controls / (float) total) < c_min_control_frac && (n_controls / total != 1)){
		std::cerr << "WARNING: Number of cases is less than " <<
				c_min_control_frac * 100 << "% of the data.  Allele frequencies"
				" for cases may be unreliable\n";
	}

	// Print a warning if rare variants will only be fixed variants
	if (1 / static_cast<float>(2 *n_controls) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in controls are rare.\n";
	}
	if (RareCaseControl && n_controls != total && 1 / static_cast<float>(2*(n_cases)) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in cases are rare.\n";
	}

}

void PopulationManager::parsePhenotypeFile(const string& filename){

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
	while(data_file.good()){
		getline(data_file, line);
		++line_no;
		trim(line);

		if(line.size() > 0 && line[0] != '#'){
			split(result, line, is_any_of(" \n\t"), boost::token_compress_on);

			if (result.size() && result.size() < 2){
				std::cerr << "WARNING: improperly formatted phenotype file "
						<< "on line " << line_no << ".\n";

			}else if(result.size()){
				if (_positions.find(result[0]) == _positions.end()){
					std::cerr << "WARNING: cannot find " << result[0] << " in VCF file.\n";
				}

				try{
					_phenotypes[result[0]] = lexical_cast<float>(result[1]);
				}catch(bad_lexical_cast&){
					_phenotypes[result[0]] = std::numeric_limits<float>::quiet_NaN();
				}
			}
		}
	}
}

float PopulationManager::getCaseAF(const Locus& loc) const{
	unordered_map<const Locus*, bitset_pair>::const_iterator itr = _genotypes.find(&loc);

	if(itr != _genotypes.end()){
		return getMAF((*itr).second, &_case_bitset);
	} else {
		return std::numeric_limits<float>::quiet_NaN();
	}
}

unsigned int PopulationManager::genotypeContribution(const Locus& loc) const{
	unordered_map<const Locus*, bitset_pair>::const_iterator itr = _genotypes.find(&loc);

	if(itr != _genotypes.end()){
		return getTotalContrib((*itr).second);
	} else {
		return 0;
	}
}

float PopulationManager::getIndivContrib(const Locus& loc, int pos, bool useWeights, const Information* const info, const Region* const reg) const{

	unordered_map<const Knowledge::Locus*, bitset_pair>::const_iterator it = _genotypes.find(&loc);

	static float weight_cache = 1;
	static const Locus* loc_cache = 0;
	float custom_weight = 1;
	unsigned int n_var = 0;

	if(it == _genotypes.end()){
		return 0;
	}

	bool g1 = (*it).second.first[pos];
	bool g2 = (*it).second.second[pos];
	n_var = (g1 && g2) ? 0 : 2*g1 + g2;

	if(c_model == DOMINANT){
		n_var = n_var > 0;
	} else if (c_model == RECESSIVE){
		n_var = n_var > 1;
	}

	// Cache the weights so we aren't wasting so much effort.
	// Also, only calculate the weight if n_var > 0 (o/w will multiply out to 0)
	if(n_var != 0 && useWeights){
		if(_use_custom_weight && info){
			custom_weight = getCustomWeight(loc, *info, reg);
		}
		if(c_use_calc_weight && loc_cache != &loc){
			loc_cache = &loc;
			weight_cache = calcWeight(loc);
		}

	}

	return n_var * weight_cache * custom_weight;
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

float PopulationManager::calcWeight(const Locus& loc) const{
	// *_a = Affected (cases), *_u = Unaffected (controls)
	// N_* = Number (population), F_* = Frequency
	int N_a, N_u, N;
	int M_a, M_u, M;
	float w_a = std::numeric_limits<float>::quiet_NaN();
	float w_u = 1;

	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotypes.find(&loc);

	if(it == _genotypes.end()){
		return 1;
	}

	dynamic_bitset<> nonmissing = ~((*it).second.first & (*it).second.second);

	N_u = (nonmissing & _control_bitset).count() ;
	N_a = (nonmissing & _case_bitset).count();
	N = N_u + N_a;

	M_u = 2*((*it).second.first & _control_bitset & nonmissing).count() +
			((*it).second.second & _control_bitset & nonmissing).count();
	M_a = 2*((*it).second.first & _case_bitset & nonmissing).count() +
			((*it).second.second & _case_bitset & nonmissing).count();
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

float PopulationManager::getCustomWeight(const Locus& loc, const Information& info, const Region* const reg) const{
	static float weight_cache = 1;
	static const Locus* loc_cache = 0;
	static const Region* reg_cache = 0;

	if(&loc != loc_cache || reg != reg_cache){
		loc_cache = &loc;
		reg_cache = reg;
		weight_cache = info.getSNPWeight(loc, reg);
	}

	return weight_cache;
}

boost::array<unsigned int, 2> PopulationManager::getBinCapacity(Bin& bin) const {

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
			capacity[0] += (nonmiss & _control_bitset).count();
			capacity[1] += (nonmiss & _case_bitset).count();
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

