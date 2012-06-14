/*
 * PopulationManager.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#include "PopulationManager.h"

#include "binmanager.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

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

const vector<bool>& PopulationManager::loadIndividuals(DataImporter& importer){
	const vector<string>& indivs = importer.getIndividualIDs();

	int size = indivs.size();
	for (int i=0; i<(size); ++i){
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
	}else if(1-(control / (float) total) < c_min_control_frac){
		std::cerr << "WARNING: Number of cases is less than " <<
				c_min_control_frac * 100 << "% of the data.  Allele frequencies"
				" for cases may be unreliable\n";
	}

	// Print a warning if rare variants will only be fixed variants
	if (1 / static_cast<float>(2 *control) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in controls are rare.\n";
	}
	if (DataImporter::RareCaseControl && control != total && 1 / static_cast<float>(2*(total - control)) > BinManager::mafCutoff){
		std::cerr << "WARNING: MAF cutoff is set so low that only variants fixed in cases are rare.\n";
	}

	return _is_control;
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




void PopulationManager::printGenotypes(ostream& os, const string& sep) const{

	unordered_map<const Locus*, bitset_pair >::const_iterator l_itr = _genotype_bitset.begin();
	unordered_map<const Locus*, bitset_pair >::const_iterator l_end = _genotype_bitset.end();

	map<string, int>::const_iterator m_itr = _positions.begin();
	map<string, int>::const_iterator m_end = _positions.end();

	map<string, float>::const_iterator pheno_status;
	map<string, float>::const_iterator pheno_end = _phenotypes.end();
	// Print the first line// TODO: format the genotype if we want to!
	os << "ID" << sep << "Status";

	while(l_itr != l_end){
		os << sep << (*l_itr).first->getID();
		++l_itr;
	}
	os << "\n";

	int pos;
	float status;
	while (m_itr != m_end){
		l_itr = _genotype_bitset.begin();
		l_end = _genotype_bitset.end();

		pos = (*m_itr).second;

		pheno_status = _phenotypes.find((*m_itr).first);
		status = -1;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second;
		}

		os << (*m_itr).first << sep << status;

		while(l_itr != l_end){
			// TODO: format the genotype if we want to!
			os << sep << (*l_itr).second.first[pos] << "/" << (*l_itr).second.second[pos];
			++l_itr;
		}

		os << "\n";
		++m_itr;
	}

}

int PopulationManager::genotypeContribution(const Locus& loc) const{
	unordered_map<const Knowledge::Locus*, int>::const_iterator itr = _genotype_sum.find(&loc);
	return (itr != _genotype_sum.end()) ? (*itr).second : 0;
}

int PopulationManager::getIndivContrib(const Locus& loc, int pos) const{
	unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator it = _genotype_bitset.find(&loc);

	if(it == _genotype_bitset.end()){
		return 0;
	}

	switch(c_model){
	case ADDITIVE:
		return (*it).second.first[pos] + (*it).second.second[pos];
	case DOMINANT:
		return (*it).second.first[pos] | (*it).second.second[pos];
	case RECESSIVE:
		return (*it).second.first[pos] & (*it).second.second[pos];
	default:
		return 0;
	}
}

array<uint, 2>& PopulationManager::getBinCapacity(Bin& bin) const{
	unordered_map<Bin*, array<uint, 2> >::const_iterator bin_itr =
			_bin_capacity.find(&bin);
	if(bin_itr == _bin_capacity.end()){
		Bin::const_locus_iterator b_itr = bin.variantBegin();
		Bin::const_locus_iterator b_end = bin.variantEnd();
		array<uint, 2> capacity;
		capacity[0] = 0;
		capacity[1] = 0;
		unordered_map<Knowledge::Locus*, array<uint, 2> >::const_iterator l_end =
				_locus_count.end();
		unordered_map<Knowledge::Locus*, array<uint, 2> >::const_iterator l_itr;
		while(b_itr != b_end){
			l_itr = _locus_count.find(*b_itr);
			if(l_itr != l_end){
				capacity[0]+=(*l_itr).second[0] / (1 + c_model != ADDITIVE);
				capacity[1]+=(*l_itr).second[1] / (1 + c_model != ADDITIVE);
			}
			++b_itr;
		}
		_bin_capacity[&bin] = capacity;
	}

	return _bin_capacity[&bin];
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
}

