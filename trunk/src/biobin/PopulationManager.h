/*
 * PopulationManager.h
 *
 *  Created on: Dec 7, 2011
 *      Author: jrw32
 */

#ifndef BIOBIN_POPULATIONMANAGER_H
#define BIOBIN_POPULATIONMANAGER_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <boost/unordered_map.hpp>

#include "knowledge/Locus.h"
#include "Bin.h"
#include "dataimporter.h"

using std::vector;
using std::string;
//using std::map;
using std::ostream;
using boost::unordered_map;

using Knowledge::Locus;

namespace BioBin{

/**
 * This is a class to manage all of the people in a given population so we can
 * view their information in a consistent way
 */
class PopulationManager{

public:
	// nothing needed for default constructor
	PopulationManager(){}

	// Loading functions
	void loadIndividuals(DataImporter& importer);
	template <class str_cont>
	void loadPhenotypes(const str_cont& phenotype_filenames);
	template <class Locus_cont>
	void loadGenotypes(const Locus_cont& dataset, DataImporter& importer);

	// Usage functions
	int genotypeContribution(const Locus& locus) const;

	// Printing functions
	template <class Bin_cont>
	void printBins(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	void printGenotypes(ostream& os, const string& sep=",") const;

private:

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	void parsePhenotypeFile(const string& filename);

	map<string, float> _phenotypes;
	map<string, int> _positions;
	map<Knowledge::Locus*, vector<short> > _genotype_map;
	unordered_map<const Knowledge::Locus*, int> _genotype_sum;

};


template <class str_cont>
void PopulationManager::loadPhenotypes(const str_cont& phenotype_files){
	vector<string>::const_iterator itr = phenotype_files.begin();
	vector<string>::const_iterator end = phenotype_files.end();

	while(itr != end){
		parsePhenotypeFile(*itr);
		++itr;
	}
}


template <class Locus_cont>
void PopulationManager::loadGenotypes(const Locus_cont& dataset, DataImporter& importer){

	typename Locus_cont::const_iterator itr = dataset.begin();
	typename Locus_cont::const_iterator end = dataset.end();

	vector<short>::const_iterator g_itr;
	vector<short>::const_iterator g_end;
	pair<uint, uint> decoded_genotype;
	int total_contrib;

	while(itr != end){
		importer.parseSNP(**itr,_genotype_map[*itr]);

		_genotype_sum[*itr] = total_contrib = 0;
		g_itr = _genotype_map[*itr].begin();
		g_end = _genotype_map[*itr].end();

		while(g_itr != g_end){
			decoded_genotype = (*itr)->decodeGenotype(*g_itr);
			total_contrib += (decoded_genotype.first != (uint)-1 && decoded_genotype.first != (*itr)->getMajorPos());
			total_contrib += (decoded_genotype.second != (uint)-1 && decoded_genotype.second != (*itr)->getMajorPos());
			++g_itr;
		}
		_genotype_sum[*itr] = total_contrib;
		++itr;
	}
}

template <class Bin_cont>
void PopulationManager::printBins(ostream& os, const Bin_cont& bins, const string& sep) const{

	typename Bin_cont::const_iterator b_itr = bins.begin();
	typename Bin_cont::const_iterator b_end = bins.end();

	// Print first line
	os << "ID" << sep << "Status";
	while(b_itr != b_end){
		os << sep << (*b_itr)->getName();
		++b_itr;
	}
	os << "\n";

	// Print second Line (totals)
	os << "Totals" << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getSize();
		++b_itr;
	}
	os << "\n";

	// Print third line (variant totals)
	os << "Variant totals" << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getVariantSize();
		++b_itr;
	}
	os << "\n";

	map<string, int>::const_iterator m_itr = _positions.begin();
	map<string, int>::const_iterator m_end = _positions.end();

	map<string, float>::const_iterator pheno_status;
	map<string, float>::const_iterator pheno_end = _phenotypes.end();

	Bin::const_locus_iterator l_itr;
	Bin::const_locus_iterator l_end;

	map<Locus*, vector<short> >::const_iterator l_pos;
	map<Locus*, vector<short> >::const_iterator l_not_found = _genotype_map.end();

	int pos;
	float status;
	int bin_count;
	pair <uint, uint> decoded_genotype;

	while (m_itr != m_end){
		b_itr = bins.begin();
		b_end = bins.end();

		pos = (*m_itr).second;

		pheno_status = _phenotypes.find((*m_itr).first);
		status = -1;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second;
		}

		os << (*m_itr).first << sep << status;

		while(b_itr != b_end){
			// Accumulate the contribution of this person in this bin
			l_itr = (*b_itr)->variantBegin();
			l_end = (*b_itr)->variantEnd();
			bin_count = 0;
			while(l_itr != l_end){
				l_pos = _genotype_map.find((*l_itr));
				if (l_pos != l_not_found){
					decoded_genotype = (*l_itr)->decodeGenotype((*l_pos).second[pos]);
					bin_count += (decoded_genotype.first != (uint)-1 && decoded_genotype.first != (*l_itr)->getMajorPos());
					bin_count += (decoded_genotype.second != (uint)-1 && decoded_genotype.second != (*l_itr)->getMajorPos());
				}
				++l_itr;
			}
			os << sep << bin_count;
			++b_itr;
		}

		os << "\n";
		++m_itr;
	}
}


}



#endif /* POPULATIONMANAGER_H_ */
