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
#include <boost/array.hpp>

#include "knowledge/Locus.h"
#include "Bin.h"
#include "dataimporter.h"

using std::vector;
using std::string;
using std::map;
using std::ostream;
using boost::array;
using boost::unordered_map;

using Knowledge::Locus;

namespace BioBin{

/**
 * This is a class to manage all of the people in a given population so we can
 * view their information in a consistent way
 */
class PopulationManager{

public:
	// A list of the available disease models.  These affect the calculation
	// of the genotype sum and an individual's contribution to a bin.
	enum DiseaseModel_ENUM { ADDITIVE, DOMINANT, RECESSIVE };

	class DiseaseModel{
	public:
		DiseaseModel() : _data(ADDITIVE){}
		DiseaseModel(const DiseaseModel_ENUM& d):_data(d){}
		operator const char*() const{
			switch(_data){
			case BioBin::PopulationManager::ADDITIVE:
				return "additive";
			case BioBin::PopulationManager::DOMINANT:
				return "dominant";
			case BioBin::PopulationManager::RECESSIVE:
				return "recessive";
			default:
				return "unknown";
			}
		}
		operator int() const{return _data;}

	private:
		DiseaseModel_ENUM _data;
	};

	// nothing needed for default constructor
	PopulationManager(){}

	// Loading functions
	const vector<bool>& loadIndividuals(DataImporter& importer);
	//template <class str_cont>
	//void loadPhenotypes(const str_cont& phenotype_filenames);
	template <class Locus_cont>
	void loadGenotypes(const Locus_cont& dataset, DataImporter& importer);

	// Usage functions
	int genotypeContribution(const Locus& locus) const;

	const vector<bool>& getControls() const {return _is_control;}

	// Printing functions
	template <class Bin_cont>
	void printBins(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	template <class Bin_cont>
	void printBinFreq(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	void printGenotypes(ostream& os, const string& sep=",") const;

	static float c_phenotype_control;
	static vector<string> c_phenotype_files;

	static float c_min_control_frac;

	static DiseaseModel c_model;

private:

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	void parsePhenotypeFile(const string& filename);

	int getIndivContrib(const Locus& loc, short genotype) const;

	array<uint, 2>& getBinCapacity(Bin& bin) const;

	map<string, float> _phenotypes;
	map<string, int> _positions;
	// _is_control can be passed to the VCF parser
	vector<bool> _is_control;
	map<Knowledge::Locus*, vector<short> > _genotype_map;
	unordered_map<const Knowledge::Locus*, int> _genotype_sum;

	unordered_map<Knowledge::Locus*, array<uint, 2> > _locus_count;
	mutable unordered_map<Bin*, array<uint, 2> > _bin_capacity;
};

template <class Locus_cont>
void PopulationManager::loadGenotypes(const Locus_cont& dataset, DataImporter& importer){

	// Get the number of people genotyped for each SNP
	importer.getNumNonMissing(dataset, _is_control, _locus_count);

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
			total_contrib += getIndivContrib(**itr, *g_itr);
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
	os << "Total Variants" << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getSize();
		++b_itr;
	}
	os << "\n";

	// Print third line (variant totals)
	os << "Total Loci" << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getVariantSize();
		++b_itr;
	}
	os << "\n";


	Bin::const_locus_iterator l_itr;
	Bin::const_locus_iterator l_end;

	unordered_map<Knowledge::Locus*, array<uint, 2> >::const_iterator loc_itr;
	unordered_map<Knowledge::Locus*, array<uint, 2> >::const_iterator loc_not_found =
			_locus_count.end();

	int locus_count = 0;
	// Print 4th + 5th lines (variant totals for cases + controls
	for(int i=0; i<2; i++){
		os << (i ? "Case" : "Control") << " Loci Totals" << sep << -1;
		b_itr = bins.begin();
		b_end = bins.end();
		while(b_itr != b_end){
			l_itr = (*b_itr)->variantBegin();
			l_end = (*b_itr)->variantEnd();
			locus_count = 0;
			while(l_itr != l_end){
				loc_itr = _locus_count.find((*l_itr));
				if (loc_itr != loc_not_found){
					locus_count += ((*loc_itr).second[i] != 0);
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
		os << (i ? "Case" : "Control") << " Bin Capacity" << sep << -1;
		b_itr = bins.begin();
		b_end = bins.end();
		while(b_itr != b_end){
			os << sep << getBinCapacity(**b_itr)[i];
			++b_itr;
		}
		os << "\n";
	}

	map<string, int>::const_iterator m_itr = _positions.begin();
	map<string, int>::const_iterator m_end = _positions.end();

	map<string, float>::const_iterator pheno_status;
	map<string, float>::const_iterator pheno_end = _phenotypes.end();

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
					//decoded_genotype = (*l_itr)->decodeGenotype((*l_pos).second[pos]);
					bin_count += getIndivContrib(**l_itr, (*l_pos).second[pos]);
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


template <class Bin_cont>
void PopulationManager::printBinFreq(ostream& os, const Bin_cont& bins, const string& sep) const{
	typename Bin_cont::const_iterator b_itr = bins.begin();
	typename Bin_cont::const_iterator b_end = bins.end();

	int n_cases = 0;
	int n_controls = 0;

	vector<bool>::const_iterator c_itr = _is_control.begin();
	vector<bool>::const_iterator c_end = _is_control.end();

	while(c_itr != c_end){
		n_controls += *c_itr;
		n_cases += !(*c_itr);
		++c_itr;
	}

	map<string, int>::const_iterator m_itr;
	map<string, int>::const_iterator m_end = _positions.end();
	int case_cont_contrib[2];

	map<Locus*, vector<short> >::const_iterator l_pos;
	map<Locus*, vector<short> >::const_iterator l_not_found = _genotype_map.end();

	os << "Bin" << sep << "Control Freq." << sep << "Case Freq.\n";

	while(b_itr != b_end){

		Bin::const_locus_iterator v_itr = (*b_itr)->variantBegin();
		Bin::const_locus_iterator v_end = (*b_itr)->variantEnd();

		case_cont_contrib[0] = 0;
		case_cont_contrib[1] = 0;
		while(v_itr != v_end){
			m_itr = _positions.begin();
			while(m_itr != m_end){
				l_pos = _genotype_map.find((*v_itr));
				if (l_pos != l_not_found){
					case_cont_contrib[!_is_control[(*m_itr).second]] +=
							getIndivContrib(**v_itr, (*l_pos).second[(*m_itr).second]);
				}
				++m_itr;
			}
			++v_itr;
		}

		os << (*b_itr)->getName() << sep;
		int capacity = 0;
		for (int i=0; i<=1; i++){
			capacity = getBinCapacity(**b_itr)[i];
			os << (capacity ? case_cont_contrib[i] / ((float) capacity) : -1);
			os << sep;
		}
		os << "\n";

		++b_itr;
	}

}

}

namespace std{
ostream& operator<<(ostream& o, const BioBin::PopulationManager::DiseaseModel& m);
}

#endif /* POPULATIONMANAGER_H_ */
