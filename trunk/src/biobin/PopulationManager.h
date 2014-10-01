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

#include "Bin.h"
#include "binmanager.h"

#include "knowledge/Locus.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/Information.h"
#include "knowledge/Region.h"

#include "vcftools/vcf_file.h"

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

	// A list of the available weighting models.  These affect the calculation
	// of the way the Madsen + Browning weight is calculated
	enum Weight_ENUM { MAX, MIN, CONTROL, OVERALL };

	class WeightModel{
	public:
		WeightModel() : _data(MAX){}
		WeightModel(const Weight_ENUM& d):_data(d){}
		operator const char*() const{
			switch(_data){
			case BioBin::PopulationManager::MAX:
				return "maximum";
			case BioBin::PopulationManager::MIN:
				return "minimum";
			case BioBin::PopulationManager::CONTROL:
				return "control";
			case BioBin::PopulationManager::OVERALL:
				return "overall";
			default:
				return "unknown";
			}
		}
		operator int() const{return _data;}

	private:
		Weight_ENUM _data;
	};

	explicit PopulationManager(const std::string& vcf_file);
	~PopulationManager(){delete vcf;}

	template <class T_cont>
	void loadLoci(T_cont& loci_out, const Knowledge::Liftover::Converter* conv);

	// Usage functions
	int genotypeContribution(const Knowledge::Locus& locus) const;

	const std::vector<bool>& getControls() const {return _is_control;}

	// Printing functions
	template <class Bin_cont>
	void printBins(std::ostream& os, const Bin_cont& bins, const Knowledge::Information& info, const std::string& sep=",") const;
	template <class Bin_cont>
	void printBinsTranspose(std::ostream& os, const Bin_cont& bins, const Knowledge::Information& info, const std::string& sep=",") const;
	template <class Bin_cont>
	void printBinFreq(std::ostream& os, const Bin_cont& bins, const std::string& sep=",") const;
	template <class Locus_cont>
	void printGenotypes(std::ostream& os, const Locus_cont& loci, const std::string& sep=",") const;

	float getCaseAF(const Knowledge::Locus& loc) const;

	static float c_phenotype_control;
	static std::vector<std::string> c_phenotype_files;

	static bool CompressedVCF;					///< gzipped file Y/N
	static bool KeepCommonLoci;
	// determine rarity of a variant by either case or control status
	static bool RareCaseControl;
	static bool OverallMajorAllele;

	static bool c_use_calc_weight;
	static bool _use_custom_weight;

	static float c_min_control_frac;

	static DiseaseModel c_model;
	static WeightModel c_weight_type;

	typedef std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<> > bitset_pair;

private:

	void printEscapedString(std::ostream& os, const std::string& toPrint, const std::string& toRepl, const std::string& replStr) const;
	std::string getEscapeString(const std::string& sep) const;

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	// Loading functions
	void loadIndividuals();
	void parsePhenotypeFile(const std::string& filename);

	float getIndivContrib(const Knowledge::Locus& loc, int position, bool useWeights = false, const Knowledge::Information* const info = NULL, const Knowledge::Region* const reg = NULL) const;
	int getTotalContrib(const Knowledge::Locus& loc) const;
	float calcBrowningWeight(unsigned long N, unsigned long M) const;
	float calcWeight(const Knowledge::Locus& loc) const;
	float getCustomWeight(const Knowledge::Locus& loc, const Knowledge::Information& info, const Knowledge::Region* const reg = NULL) const;

	float getMAF(const std::vector<int>& allele_count, unsigned int nmcc) const;

	boost::array<unsigned int, 2> getBinCapacity(Bin& bin) const;

	std::map<std::string, float> _phenotypes;
	std::map<std::string, int> _positions;

	// _is_control can be passed to the VCF parser
	std::vector<bool> _is_control;
	boost::dynamic_bitset<> _control_bitset;
	int n_controls;

	boost::unordered_map<const Knowledge::Locus*, bitset_pair > _genotype_bitset;

	// A map to keep track of where in the file a locus resides
	boost::unordered_map<Knowledge::Locus*, int> _locus_position;

	boost::unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> > _locus_count;

	VCF::vcf_file* vcf;
};

template <class Bin_cont>
void PopulationManager::printBinsTranspose(std::ostream& os, const Bin_cont& bins, const Knowledge::Information& info, const std::string& sep) const{

	_use_custom_weight = (info.c_weight_files.size() > 0);

	std::string sep_repl = getEscapeString(sep);

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

	std::map<std::string, int>::const_iterator m_itr = _positions.begin();
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
	std::map<std::string, float>::const_iterator pheno_status;
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
	boost::unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator loc_itr;

	while(b_itr != bins.end()){
		// Print bin name
		printEscapedString(os, (*b_itr)->getName(), sep, sep_repl);

		// print total var
		os << sep << (*b_itr)->getSize();

		// print total loci
		os << sep << (*b_itr)->getVariantSize();

		// print case/control loci
		boost::array<unsigned int, 2> num_loci;
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
		boost::array<unsigned int, 2> capacity = getBinCapacity(**b_itr);
		os << sep << capacity[0] << sep << capacity[1];

		// print for each person
		float bin_count = 0;
		m_itr = _positions.begin();
		while (m_itr != _positions.end()) {
			l_itr = (*b_itr)->variantBegin();
			bin_count = 0;
			while (l_itr != (*b_itr)->variantEnd()) {
				bin_count += getIndivContrib(**l_itr, (*m_itr).second, true, &info, (*b_itr)->getRegion());
				++l_itr;
			}
			os << sep << std::setprecision(4) << bin_count;
			++m_itr;
		}
		os << "\n";
		++b_itr;
	}
}

template <class Bin_cont>
void PopulationManager::printBins(std::ostream& os, const Bin_cont& bins, const Knowledge::Information& info, const std::string& sep) const{
	std::string sep_repl = getEscapeString(sep);

	_use_custom_weight = (info.c_weight_files.size() > 0);

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

	boost::unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator loc_itr;
	boost::unordered_map<const Knowledge::Locus*, boost::array<unsigned short, 2> >::const_iterator loc_not_found =
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

	std::map<std::string, int>::const_iterator m_itr = _positions.begin();
	std::map<std::string, int>::const_iterator m_end = _positions.end();

	std::map<std::string, float>::const_iterator pheno_status;
	std::map<std::string, float>::const_iterator pheno_end = _phenotypes.end();

	int pos;
	float status;
	float bin_count;

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
				bin_count += getIndivContrib(**l_itr, pos, true, &info, (*b_itr)->getRegion());
				++l_itr;
			}
			os << sep << std::setprecision(4) << bin_count;
			++b_itr;
		}

		os << "\n";
		++m_itr;
	}
}


template <class Bin_cont>
void PopulationManager::printBinFreq(std::ostream& os, const Bin_cont& bins, const std::string& sep) const{

	std::string sep_repl = getEscapeString(sep);

	typename Bin_cont::const_iterator b_itr = bins.begin();
	typename Bin_cont::const_iterator b_end = bins.end();

	std::map<std::string, int>::const_iterator m_itr;
	std::map<std::string, int>::const_iterator m_end = _positions.end();
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
						getIndivContrib(**v_itr, (*m_itr).second, false);
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
void PopulationManager::printGenotypes(std::ostream& os, const Locus_cont& loci, const std::string& sep) const{

	string sep_repl = getEscapeString(sep);

	typename Locus_cont::const_iterator l_itr = loci.begin();

	std::map<std::string, int>::const_iterator m_itr = _positions.begin();
	std::map<std::string, int>::const_iterator m_end = _positions.end();

	std::map<std::string, float>::const_iterator pheno_status;
	std::map<std::string, float>::const_iterator pheno_end = _phenotypes.end();

	boost::unordered_map<const Knowledge::Locus*, bitset_pair>::const_iterator g_itr;
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
	std::set<std::string> unknownChromosomes;
	unsigned int totalSiteCount	= vcf->N_entries;

	// Predefine everything so that the loop below can be as tight as possible
	std::vector<std::pair<int, int> > genotype_pairs;
	genotype_pairs.resize(size);
	bitset_pair gen_bits(std::make_pair(boost::dynamic_bitset<>(size),boost::dynamic_bitset<>(size)));
	std::pair<int,int> genotype;
	boost::array<unsigned short, 2> nm;
	unsigned int alleleCount = 0;
	std::string line;
	boost::array<std::vector<int>, 2> alleleCounts;
	std::pair<boost::unordered_map<const Knowledge::Locus*, bitset_pair>::iterator, bool> gen_pair;
	boost::unordered_map<const Knowledge::Locus*, bitset_pair>::iterator gen_itr;
	VCF::vcf_entry entry(vcf->N_indv);

	for (unsigned int i=0; i<totalSiteCount; i++) {
		vcf->get_vcf_entry(i, line);
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

		for (unsigned int j=0; j<vcf->N_indv; j++) {
			entry.get_indv_GENOTYPE_ids(j, genotype);

			if(genotype.first != -1){
				++alleleCounts[!_control_bitset[j]][genotype.first];
			}
			if(genotype.second != -1){
				++alleleCounts[!_control_bitset[j]][genotype.second];
			}

			if((genotype.first==-1)^(genotype.second==-1)){
				std::cerr << "WARNING: Genotype partially missing at chromosome "
						<< entry.get_CHROM() <<", position " << entry.get_POS()
						<< std::endl;
			}

			nm[!_control_bitset[j]] += (genotype.first != -1);
			nm[!_control_bitset[j]] += (genotype.second != -1);

			// I need to save this to determine if it is in fact "minor"
			genotype_pairs[j] = genotype;

		}

		// To deal with the default parameter, we really just want to use
		// vcf->include_indivs, but NOOOO C++ has to be a pain

		float controlMAF = getMAF(alleleCounts[0], nm[0]);
		float caseMAF = getMAF(alleleCounts[1], nm[1]);

		bool is_rare = ( controlMAF <= BinManager::mafCutoff && controlMAF >= BinManager::mafThreshold) ||
				(RareCaseControl && nm[1] >= 0 && caseMAF <= BinManager::mafCutoff && caseMAF >= BinManager::mafThreshold);

		// If the minimum bin size is > 0 and this locus has no minor alleles,
		// just keep going!

		// Note that if the MAF is < 0, we have no data for that population, so
		// it could contribute nothing to the analysis
		bool drop_loci = controlMAF <= 0 && caseMAF <= 0 && BinManager::MinBinSize > 0;

		if(KeepCommonLoci || (is_rare && !drop_loci) ){
			Knowledge::Locus* loc = new Knowledge::Locus(entry.get_CHROM(),entry.get_POS(),is_rare,entry.get_ID());

			if(conv){
				Knowledge::Locus* new_loc = conv->convertLocus(*loc);
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
				string currMajor = "";
				float currMax = 0;
				float currFreq = 0;

				loci_out.insert(loci_out.end(), loc);
				_locus_position.insert(make_pair(loc, i));

				// Add ability to order the alleles by overall population here!
				float freq = (nm[0] != 0 ? (alleleCounts[0][0] / static_cast<float>(nm[0])) : 0);
				int nm_sum = nm[0] + nm[1];
				if (OverallMajorAllele){
					currMax = (nm_sum != 0 ? (alleleCounts[0][0] + alleleCounts[1][0]) / static_cast<float>(nm_sum) : 0);
					currMajor = entry.get_REF();
				}


				loc->addAllele(entry.get_REF(), freq);
				for (uint n = 1; n<alleleCount; n++){
					freq = (nm[0] != 0 ? alleleCounts[0][n] / static_cast<float>(nm[0]) : -1);
					if (OverallMajorAllele){
						currFreq = (nm_sum != 0 ? (alleleCounts[0][n] + alleleCounts[1][n]) / static_cast<float>(nm_sum) : -1);
						if(currFreq > currMax){
							currMax = currFreq;
							currMajor = entry.get_ALT_allele(n-1);
						}
					}
					loc->addAllele(entry.get_ALT_allele(n-1), freq);
				}
				if(OverallMajorAllele){
					loc->setMajorAllele(currMajor);
				}

				_locus_count.insert(make_pair(loc, nm));
				gen_pair = _genotype_bitset.insert(
						std::make_pair(loc,
								std::make_pair(boost::dynamic_bitset<>(size),
										boost::dynamic_bitset<>(size))));
				gen_itr = gen_pair.first;

				for (unsigned int j=0; j<vcf->N_indv; ++j) {

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
istream& operator>>(istream& in, BioBin::PopulationManager::WeightModel& model_out);
ostream& operator<<(ostream& o, const BioBin::PopulationManager::WeightModel& m);
}

#endif /* POPULATIONMANAGER_H_ */
