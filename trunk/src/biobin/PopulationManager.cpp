/*
 * PopulationManager.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#include "PopulationManager.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using boost::algorithm::split;
using boost::algorithm::is_any_of;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using Knowledge::Locus;

namespace BioBin{

void PopulationManager::loadIndividuals(DataImporter& importer){
	const vector<string>& indivs = importer.getIndividualIDs();

	int size = indivs.size();
	for (int i=-1; i<(size-1); ++i){
		_positions[indivs[i]] = i;
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
	while(data_file.good()){
		getline(data_file, line);
		split(result, line, is_any_of(" \n\t"));

		if (result.size() < 2){
			std::cerr << "WARNING: improperly formatted phenotype file";
		}

		try{
			_phenotypes[result[0]] = lexical_cast<int>(result[1]);
		}catch(bad_lexical_cast&){
			_phenotypes[result[0]] = -1;
		}
	}
}




void PopulationManager::printGenotypes(ostream& os, const string& sep) const{

	map<Locus*, vector<short> >::const_iterator l_itr = _genotype_map.begin();
	map<Locus*, vector<short> >::const_iterator l_end = _genotype_map.end();

	map<string, int>::const_iterator m_itr = _positions.begin();
	map<string, int>::const_iterator m_end = _positions.end();

	map<string, int>::const_iterator pheno_status;
	map<string, int>::const_iterator pheno_end = _phenotypes.end();
	// Print the first line// TODO: format the genotype if we want to!
	os << "ID" << sep << "Status";

	while(l_itr != l_end){
		os << sep << (*l_itr).first->getID();
		++l_itr;
	}
	os << "\n";

	int pos;
	int status;
	while (m_itr != m_end){
		l_itr = _genotype_map.begin();
		l_end = _genotype_map.end();

		pos = (*m_itr).second;

		pheno_status = _phenotypes.find((*m_itr).first);
		status = -1;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second;
		}

		os << (*m_itr).first << sep << status;

		while(l_itr != l_end){
			// TODO: format the genotype if we want to!
			os << sep << (*l_itr).second[pos];
			++l_itr;
		}

		os << "\n";
		++m_itr;
	}

}





}


