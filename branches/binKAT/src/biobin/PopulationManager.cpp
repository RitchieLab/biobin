/*
 * PopulationManager.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#include "PopulationManager.h"

#include "binapplication.h"
#include "binmanager.h"

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
string PopulationManager::c_phenotype_file = "";
float PopulationManager::c_min_control_frac = 0.125;
PopulationManager::DiseaseModel PopulationManager::c_model =
		PopulationManager::ADDITIVE;
PopulationManager::WeightModel PopulationManager::c_weight_type =
		PopulationManager::MAX;

bool PopulationManager::RareCaseControl = true;
bool PopulationManager::c_use_calc_weight = false;
bool PopulationManager::_use_custom_weight = false;
bool PopulationManager::NoSummary = false;

PopulationManager::PopulationManager(const string& vcf_fn) : _vcf_fn(vcf_fn){

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

	bool allControl = false;
	if(c_phenotype_file != ""){
		parsePhenotypeFile(c_phenotype_file);
	} else {
		allControl = true;
	}

	// Now, we go through and determine who is a case and who is a control
	boost::unordered_map<string, int>::const_iterator p_itr = _positions.begin();
	boost::unordered_map<string, int>::const_iterator p_end = _positions.end();
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
	bool header_read = false;
	while(data_file.good()){
		getline(data_file, line);
		++line_no;
		trim(line);

		if(line.size() > 0){

			if(!header_read && line[0] == '#'){
				split(result, line, is_any_of(" \n\t"), boost::token_compress_on);
				for(unsigned int i=1; i<result.size(); i++){
					_pheno_names.push_back(result[i]);
				}
			} else if(line[0] != '#'){
				split(result, line, is_any_of(" \n\t"), boost::token_compress_on);

				// if we haven't read the header, just assign arbitrary phenotype names
				if(!header_read){
					if(result.size() > 2){
						std::cerr << "WARNING: No header given for multiple phenotypes, assigning sequential phenotype names" << std::endl;
						for(unsigned int i=1; i<result.size(); i++){
							_pheno_names.push_back("pheno_" + boost::lexical_cast<string>(i));
						}
					}else{
						_pheno_names.push_back("");
					}
				}

				if (result.size() != _pheno_names.size() + 1){
					std::cerr << "WARNING: improperly formatted phenotype file "
							  << "on line " << line_no << ", ignoring." << std::endl;
				}else if (_positions.find(result[0]) == _positions.end()){
						std::cerr << "WARNING: cannot find " << result[0] << " in VCF file." << std::endl;
				} else {
					_phenos[result[0]].reserve(_pheno_names.size());
					for(unsigned int i=1; i<result.size(); i++){
						try{
							_phenos[result[0]].push_back(lexical_cast<float>(result[i]));
						}catch(bad_lexical_cast&){
							_phenos[result[0]].push_back(std::numeric_limits<float>::quiet_NaN());
						}
					}
				}
			}
			header_read = true;
		}
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

float PopulationManager::getIndivContrib(const Locus& loc, int pos, const bitset_pair& status, bool useWeights, const Information* const info, const Region* const reg) const{

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
			weight_cache = calcWeight(loc, status);
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

