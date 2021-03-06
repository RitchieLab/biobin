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

bool PopulationManager::CompressedVCF = false;
bool PopulationManager::KeepCommonLoci = true;
bool PopulationManager::RareCaseControl = true;
bool PopulationManager::OverallMajorAllele = true;
bool PopulationManager::c_use_calc_weight = false;
bool PopulationManager::_use_custom_weight = false;

PopulationManager::PopulationManager(const string& vcf_fn) : vcf(NULL){
	if(BinApplication::s_run_normal){
		vcf = new VCF::vcf_file(vcf_fn, CompressedVCF);
	}
}

void PopulationManager::loadIndividuals(){
	const vector<string>& indivs = vcf->indv;
	//importer.getIndividualIDs();

	int size = indivs.size();
	for (int i=0; i<(size); i++){
		_positions[indivs[i]] = i;
	}

	// By default, everyone is a control who is not found in a phenotype file
	_is_control = vector<bool>(size, true);
	_control_bitset = dynamic_bitset<>(size);
	_control_bitset.set();

	// OK, now iterate through the phenotype files and load them up
	vector<string>::const_iterator itr = c_phenotype_files.begin();
	vector<string>::const_iterator end = c_phenotype_files.end();

	while(itr != end){
		parsePhenotypeFile(*itr);
		++itr;
	}

	// Now, we go through and determine who is a case and who is a control
	map<string, int>::const_iterator p_itr = _positions.begin();
	map<string, int>::const_iterator p_end = _positions.end();
	map<string, float>::const_iterator pheno_itr;
	map<string, float>::const_iterator pheno_not_found = _phenotypes.end();
	int total = 0;
	int control = 0;
	while(p_itr != p_end){
		pheno_itr = _phenotypes.find((*p_itr).first);
		++total;
		if (pheno_itr == pheno_not_found || (*pheno_itr).second == c_phenotype_control){
			++control;
		}else{
			int pos = (*p_itr).second;
			_is_control[pos] = false;
			_control_bitset.reset(pos);
		}
		++p_itr;
	}

	// If we don't have enough controls, print a warning and make everyone a control
	if ((control / (float) total) < c_min_control_frac){
		std::cerr << "WARNING: Number of controls is less than " <<
				c_min_control_frac * 100 << "% of the data.  Using all individuals as controls\n";
		_is_control = vector<bool>(size, true);
		_control_bitset.set();
		control = total;
	}else if(1-(control / (float) total) < c_min_control_frac && (control / total != 1)){
		std::cerr << "WARNING: Number of cases is less than " <<
				c_min_control_frac * 100 << "% of the data.  Allele frequencies"
				" for cases may be unreliable\n";
	}

	// Print a warning if rare variants will only be fixed variants
	if (1 / static_cast<float>(2 *control) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in controls are rare.\n";
	}
	if (RareCaseControl && control != total && 1 / static_cast<float>(2*(total - control)) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in cases are rare.\n";
	}

	std::vector<bool>::const_iterator c_itr = _is_control.begin();
	n_controls = 0;
	while(c_itr != _is_control.end()){
		n_controls += *c_itr;
		++c_itr;
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
					_phenotypes[result[0]] = -1;
				}
			}
		}
	}
}

float PopulationManager::getCaseAF(const Locus& loc) const{
	float ret_val = -1;

	unordered_map<const Locus*, boost::array<unsigned short, 2> >::const_iterator count_itr = _locus_count.find(&loc);
	if (count_itr != _locus_count.end() && (*count_itr).second[1] != 0){
		unordered_map<const Locus*, bitset_pair>::const_iterator bitset_itr = _genotype_bitset.find(&loc);
		if(bitset_itr != _genotype_bitset.end()){
			ret_val =
					(((*bitset_itr).second.first & (~_control_bitset)).count()
							+ ((*bitset_itr).second.second & (~_control_bitset)).count())
							/ static_cast<float>((*count_itr).second[1]);
		}
	}

	return ret_val;
}

int PopulationManager::genotypeContribution(const Locus& loc) const{

	return getTotalContrib(loc);

}

float PopulationManager::getIndivContrib(const Locus& loc, int pos, bool useWeights, const Information* const info, const Region* const reg) const{
	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotype_bitset.find(&loc);

	int n_var = 0;
	static float weight_cache = 1;
	static const Locus* loc_cache = 0;
	float custom_weight = 1;

	if(it == _genotype_bitset.end()){
		return 0;
	}

	switch(c_model){
	case ADDITIVE:
		n_var = (*it).second.first[pos] + (*it).second.second[pos];
		break;
	case DOMINANT:
		n_var = (*it).second.first[pos] | (*it).second.second[pos];
		break;
	case RECESSIVE:
		n_var = (*it).second.first[pos] & (*it).second.second[pos];
		break;
	default:
		return 0;
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

int PopulationManager::getTotalContrib(const Locus& loc) const{
	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotype_bitset.find(&loc);

	int n_var = 0;

	if(it == _genotype_bitset.end()){
		return 0;
	}

	switch(c_model){
	case ADDITIVE:
		n_var = (*it).second.first.count() + (*it).second.second.count();
		break;
	case DOMINANT:
		n_var = ((*it).second.first | (*it).second.second).count();
		break;
	case RECESSIVE:
		n_var = ((*it).second.first & (*it).second.second).count();
		break;
	default:
		return 0;
	}

	return n_var;
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

	unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator l_itr = _locus_count.find(&loc);
	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotype_bitset.find(&loc);

	if(l_itr == _locus_count.end() || it == _genotype_bitset.end()){
		return 1;
	}

	N_u = (*l_itr).second[0]/2;
	N_a = (*l_itr).second[1]/2;
	N = N_u + N_a;

	M_u = ((*it).second.first & _control_bitset).count() +
			((*it).second.second & _control_bitset).count();
	M_a = ((*it).second.first & (~_control_bitset)).count() +
			((*it).second.second & (~_control_bitset)).count();
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
	unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator l_end =
			_locus_count.end();
	unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator l_itr;
	while (b_itr != b_end) {

		l_itr = _locus_count.find(*b_itr);
		if (l_itr != l_end) {
			capacity[0] += (*l_itr).second[0];
			capacity[1] += (*l_itr).second[1];
		}
		++b_itr;
	}

	capacity[0] /= (1 + c_model != ADDITIVE);
	capacity[1] /= (1 + c_model != ADDITIVE);

	return capacity;
}

float PopulationManager::getMAF(const vector<int>& allele_count, uint nmcc) const{
	if(nmcc == 0){
		return -1;
	}else{
		vector<int> ordered_allele_count(allele_count.begin(), allele_count.end());
		sort(ordered_allele_count.begin(), ordered_allele_count.end());
		float ret_val = 0;
		if (ordered_allele_count.size() > 1){
			ret_val = (*(++ordered_allele_count.rbegin())) / ((float) nmcc);
		}
		return ret_val;
	}
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

