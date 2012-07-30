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
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

#include "knowledge/Locus.h"
#include "knowledge/liftover/Converter.h"
#include "Bin.h"
#include "binmanager.h"
//#include "dataimporter.h"
#include "vcftools/vcf_file.h"

using std::fill;
using std::vector;
using std::string;
using std::map;
using std::ostream;
using boost::array;
using boost::unordered_map;
using boost::dynamic_bitset;

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
	//PopulationManager(){}

	explicit PopulationManager(const string& vcf_file);



	//template <class Locus_cont>
	//void loadGenotypes(const Locus_cont& dataset, DataImporter& importer);

	template <class T_cont>
	void loadLoci(T_cont& loci_out, const Knowledge::Liftover::Converter* conv);

	// Usage functions
	int genotypeContribution(const Locus& locus) const;

	const vector<bool>& getControls() const {return _is_control;}

	// Printing functions
	template <class Bin_cont>
	void printBins(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	template <class Bin_cont>
	void printBinsTranspose(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	template <class Bin_cont>
	void printBinFreq(ostream& os, const Bin_cont& bins, const string& sep=",") const;
	template <class Locus_cont>
	void printGenotypes(ostream& os, const Locus_cont& loci, const string& sep=",") const;

	float getCaseAF(const Locus& loc) const;

	static float c_phenotype_control;
	static vector<string> c_phenotype_files;

	static bool CompressedVCF;					///< gzipped file Y/N
	static bool KeepCommonLoci;
	// determine rarity of a variant by either case or control status
	static bool RareCaseControl;
	static bool OverallMajorAllele;

	static float c_min_control_frac;

	static DiseaseModel c_model;

	typedef pair<dynamic_bitset<>, dynamic_bitset<> > bitset_pair;

private:

	void printEscapedString(ostream& os, const string& toPrint, const string& toRepl, const string& replStr) const;
	string getEscapeString(const string& sep) const;

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	// Loading functions
	void loadIndividuals();
	void parsePhenotypeFile(const string& filename);

	int getIndivContrib(const Locus& loc, int position) const;
	int getTotalContrib(const Locus& loc) const;

	float getMAF(const vector<int>& allele_count, uint nmcc) const;

	array<unsigned int, 2> getBinCapacity(Bin& bin) const;

	map<string, float> _phenotypes;
	map<string, int> _positions;

	// _is_control can be passed to the VCF parser
	vector<bool> _is_control;
	dynamic_bitset<> _control_bitset;

	unordered_map<const Knowledge::Locus*, bitset_pair > _genotype_bitset;

	// A map to keep track of where in the file a locus resides
	unordered_map<Knowledge::Locus*, int> _locus_position;

	unordered_map<const Knowledge::Locus*, array<unsigned short, 2> > _locus_count;

	mutable VCF::vcf_file vcf;
};

template <class Bin_cont>
void PopulationManager::printBinsTranspose(ostream& os, const Bin_cont& bins, const string& sep) const{
	string sep_repl = getEscapeString(sep);

	// Print 1st line
	printEscapedString(os, "Bin Name", sep, sep_repl);
	os << sep;
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

	map<string, int>::const_iterator m_itr = _positions.begin();
	while(m_itr != _positions.end()){
		os << sep;
		printEscapedString(os, (*m_itr).first, sep, sep_repl);
		++m_itr;
	}

	os << "\n";

	printEscapedString(os, "Status", sep, sep_repl);
	os << sep <<
			-1 << sep << -1 << sep << -1 << sep << -1 << sep << -1 << sep << -1;

	m_itr = _positions.begin();
	map<string, float>::const_iterator pheno_status;
	while(m_itr != _positions.end()){
		float status = -1;
		pheno_status = _phenotypes.find((*m_itr).first);
		if (pheno_status != _phenotypes.end()){
			status = (*pheno_status).second;
		}
		os << sep << status;
		++m_itr;
	}

	os << "\n";

	typename Bin_cont::const_iterator b_itr = bins.begin();
	Bin::const_locus_iterator l_itr;
	unordered_map<const Knowledge::Locus*, array<unsigned short, 2> >::const_iterator loc_itr;

	while(b_itr != bins.end()){
		// Print bin name
		printEscapedString(os, (*b_itr)->getName(), sep, sep_repl);

		// print total var
		os << sep << (*b_itr)->getSize();

		// print total loci
		os << sep << (*b_itr)->getVariantSize();

		// print case/control loci
		array<unsigned int, 2> num_loci;
		num_loci[0] = 0;
		num_loci[1] = 0;
		l_itr = (*b_itr)->variantBegin();
		while(l_itr != (*b_itr)->variantEnd()){
			loc_itr = _locus_count.find((*l_itr));
			if (loc_itr != _locus_count.end()){
				num_loci[0] += (*loc_itr).second[0] != 0;
				num_loci[1] += (*loc_itr).second[1] != 0;
			}
			++l_itr;
		}

		os << sep << num_loci[0] << sep << num_loci[1];

		// print case/control capacity
		array<unsigned int, 2> capacity = getBinCapacity(**b_itr);
		os << sep << capacity[0] << sep << capacity[1];

		// print for each person
		int bin_count = 0;
		m_itr = _positions.begin();
		while (m_itr != _positions.end()) {
			l_itr = (*b_itr)->variantBegin();
			bin_count = 0;
			while (l_itr != (*b_itr)->variantEnd()) {
				bin_count += getIndivContrib(**l_itr, (*m_itr).second);
				++l_itr;
			}
			os << sep << bin_count;
			++m_itr;
		}
		os << "\n";
		++b_itr;
	}
}

template <class Bin_cont>
void PopulationManager::printBins(ostream& os, const Bin_cont& bins, const string& sep) const{
	string sep_repl = getEscapeString(sep);

	typename Bin_cont::const_iterator b_itr = bins.begin();
	typename Bin_cont::const_iterator b_end = bins.end();

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

	// Print second Line (totals)
	printEscapedString(os, "Total Variants", sep, sep_repl);
	os << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getSize();
		++b_itr;
	}
	os << "\n";

	// Print third line (variant totals)
	printEscapedString(os, "Total Loci", sep, sep_repl);
	os << sep << -1;
	b_itr = bins.begin();
	b_end = bins.end();
	while(b_itr != b_end){
		os << sep << (*b_itr)->getVariantSize();
		++b_itr;
	}
	os << "\n";


	Bin::const_locus_iterator l_itr;
	Bin::const_locus_iterator l_end;

	unordered_map<const Knowledge::Locus*, array<unsigned short, 2> >::const_iterator loc_itr;
	unordered_map<const Knowledge::Locus*, array<unsigned short, 2> >::const_iterator loc_not_found =
			_locus_count.end();

	int locus_count = 0;
	// Print 4th + 5th lines (variant totals for cases + controls
	for(int i=0; i<2; i++){
		printEscapedString(os, string(i ? "Case" : "Control") + " Loci Totals", sep, sep_repl);
		os << sep << -1;
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
		printEscapedString(os, string(i ? "Case" : "Control") + " Bin Capacity", sep, sep_repl);
		os << sep << -1;
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

	int pos;
	float status;
	int bin_count;

	while (m_itr != m_end){
		b_itr = bins.begin();
		b_end = bins.end();

		pos = (*m_itr).second;
		pheno_status = _phenotypes.find((*m_itr).first);
		status = -1;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second;
		}

		printEscapedString(os, (*m_itr).first, sep, sep_repl);

		os << sep << status;

		while(b_itr != b_end){
			// Accumulate the contribution of this person in this bin
			l_itr = (*b_itr)->variantBegin();
			l_end = (*b_itr)->variantEnd();
			bin_count = 0;
			while(l_itr != l_end){
				bin_count += getIndivContrib(**l_itr, pos);
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
	string sep_repl = getEscapeString(sep);

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

	printEscapedString(os, "Bin", sep, sep_repl);
	os << sep;
	printEscapedString(os, "Control Freq.", sep, sep_repl);
	os << sep;
	printEscapedString(os, "Case Freq.", sep, sep_repl);
	os << "\n";

	while(b_itr != b_end){

		Bin::const_locus_iterator v_itr = (*b_itr)->variantBegin();
		Bin::const_locus_iterator v_end = (*b_itr)->variantEnd();

		case_cont_contrib[0] = 0;
		case_cont_contrib[1] = 0;
		while(v_itr != v_end){
			m_itr = _positions.begin();
			while(m_itr != m_end){
				case_cont_contrib[!_is_control[(*m_itr).second]] +=
						getIndivContrib(**v_itr, (*m_itr).second);
				++m_itr;
			}
			++v_itr;
		}

		printEscapedString(os, (*b_itr)->getName(), sep, sep_repl);
		os << sep;

		short capacity = 0;
		for (int i=0; i<=1; i++){
			capacity = getBinCapacity(**b_itr)[i];
			os << (capacity ? case_cont_contrib[i] / ((float) capacity) : -1);
			os << sep;
		}
		os << "\n";

		++b_itr;
	}

}

template <class Locus_cont>
void PopulationManager::printGenotypes(ostream& os, const Locus_cont& loci, const string& sep) const{

	string sep_repl = getEscapeString(sep);

	typename Locus_cont::const_iterator l_itr = loci.begin();

	map<string, int>::const_iterator m_itr = _positions.begin();
	map<string, int>::const_iterator m_end = _positions.end();

	map<string, float>::const_iterator pheno_status;
	map<string, float>::const_iterator pheno_end = _phenotypes.end();

	unordered_map<const Locus*, bitset_pair>::const_iterator g_itr;
	// Print the first line// TODO: format the genotype if we want to!
	printEscapedString(os, "ID", sep, sep_repl);
	os << sep;
	printEscapedString(os, "Status", sep, sep_repl);

	while(l_itr != loci.end()){
		os << sep;
		printEscapedString(os, (*l_itr)->getID(), sep, sep_repl);
		++l_itr;
	}
	os << "\n";

	int pos;
	float status;
	while (m_itr != m_end){
		l_itr = loci.begin();

		pos = (*m_itr).second;

		pheno_status = _phenotypes.find((*m_itr).first);
		status = -1;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second;
		}

		printEscapedString(os, (*m_itr).first, sep, sep_repl);
		os << sep << status;

		while(l_itr != loci.end()){
			g_itr = _genotype_bitset.find(*l_itr);
			if(g_itr != _genotype_bitset.end()){
				// TODO: format the genotype if we want to!
				os << sep << (*g_itr).second.first[pos] << "/" << (*g_itr).second.second[pos];
			}else{
				os << sep << "?/?";
			}
			++l_itr;
		}

		os << "\n";
		++m_itr;
	}

}

template<class T_cont>
void PopulationManager::loadLoci(T_cont& loci_out, const Knowledge::Liftover::Converter* conv=0){

	loadIndividuals();
	int size = _is_control.size();
	set<string> unknownChromosomes;
	uint totalSiteCount	= vcf.N_entries;

	// Predefine everything so that the loop below can be as tight as possible
	vector<pair<int, int> > genotype_pairs;
	genotype_pairs.resize(size);
	bitset_pair gen_bits(make_pair(dynamic_bitset<>(size),dynamic_bitset<>(size)));
	pair<int,int> genotype;
	array<unsigned short, 2> nm;
	uint alleleCount = 0;
	string line;
	array<vector<int>, 2> alleleCounts;
	pair<unordered_map<const Locus*, bitset_pair>::iterator, bool> gen_pair;
	unordered_map<const Locus*, bitset_pair>::iterator gen_itr;
	VCF::vcf_entry entry(vcf.N_indv);

	for (uint i=0; i<totalSiteCount; i++) {
		vcf.get_vcf_entry(i, line);
		entry.reset(line);
		entry.parse_basic_entry(true);

		entry.parse_genotype_entries(true, false, false, false);
		alleleCount = entry.get_N_alleles();

		nm[0] = 0;
		nm[1] = 0;
		alleleCounts[0].resize(alleleCount);
		fill(alleleCounts[0].begin(), alleleCounts[0].end(), 0);
		alleleCounts[1].resize(alleleCount);
		fill(alleleCounts[1].begin(), alleleCounts[1].end(), 0);

		for (uint j=0; j<vcf.N_indv; j++) {
			entry.get_indv_GENOTYPE_ids(j, genotype);

			if(genotype.first != -1){
				++alleleCounts[!_control_bitset[j]][genotype.first];
			}
			if(genotype.second != -1){
				++alleleCounts[!_control_bitset[j]][genotype.second];
			}

			nm[!_control_bitset[j]] += (genotype.first != -1);
			nm[!_control_bitset[j]] += (genotype.second != -1);

			// I need to save this to determine if it is in fact "minor"
			genotype_pairs[j] = genotype;

		}

		// To deal with the default parameter, we really just want to use
		// vcf.include_indivs, but NOOOO C++ has to be a pain

		bool is_rare = (getMAF(alleleCounts[0], nm[0]) < BinManager::mafCutoff) ||
				(RareCaseControl && nm[1] > 0 && getMAF(alleleCounts[1], nm[1]) < BinManager::mafCutoff);

		if(KeepCommonLoci || is_rare ){
			Locus* loc = new Locus(entry.get_CHROM(),entry.get_POS(),is_rare,entry.get_ID());

			if(conv){
				Locus* new_loc = conv->convertLocus(*loc);
				if (! new_loc){
					// Add to the
					delete loc;
					loc = 0;
				}else{
					delete loc;
					loc = new_loc;
				}
			}

			if(loc){
				loci_out.insert(loci_out.end(), loc);
				_locus_position.insert(make_pair(loc, i));

				// Add ability to order the alleles by overall population here!
				float freq = (nm[0] != 0 ? (alleleCounts[0][0] / static_cast<float>(nm[0])) : 0);
				int nm_sum = nm[0] + nm[1];
				if (OverallMajorAllele){
					freq = (nm_sum != 0 ? (alleleCounts[0][0] + alleleCounts[1][0]) / static_cast<float>(nm_sum) : 0);
				}


				loc->addAllele(entry.get_REF(), freq);
				for (uint n = 1; n<alleleCount; n++){
					freq = (nm[0] != 0 ? alleleCounts[0][n] / static_cast<float>(nm[0]) : -1);
					if (OverallMajorAllele){
						freq = (nm_sum != 0 ? (alleleCounts[0][n] + alleleCounts[1][n]) / static_cast<float>(nm_sum) : -1);
					}
					loc->addAllele(entry.get_ALT_allele(n-1), freq);
				}

				_locus_count.insert(make_pair(loc, nm));
				gen_pair = _genotype_bitset.insert(
						make_pair(loc,
								make_pair(dynamic_bitset<>(size),
										dynamic_bitset<>(size))));
				gen_itr = gen_pair.first;

				for (uint j=0; j<vcf.N_indv; ++j) {

					//MOVE THIS!!
					(*gen_itr).second.first[j] = genotype_pairs[j].first != -1 && static_cast<unsigned short>(genotype_pairs[j].first) != loc->getMajorPos();
					(*gen_itr).second.second[j] = genotype_pairs[j].second != -1 && static_cast<unsigned short>(genotype_pairs[j].second) != loc->getMajorPos();
				}
			}
		}
	}
}

}

namespace std{
istream& operator>>(istream& in, BioBin::PopulationManager::DiseaseModel& model_out);
ostream& operator<<(ostream& o, const BioBin::PopulationManager::DiseaseModel& m);
}

#endif /* POPULATIONMANAGER_H_ */
