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
#include <iomanip>
#include <algorithm>
#include <stdexcept>

#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "Bin.h"

#include "util/ICompressedFile.h"
#include "util/string_ref.hpp"

#include "knowledge/Locus.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/Information.h"
#include "knowledge/Region.h"

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

	typedef std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<> > bitset_pair;

	class Phenotype{
	public:
		Phenotype(unsigned int idx, const bitset_pair& status) : _idx(idx), _status(&status) {}
		const bitset_pair& getStatus() const {return *_status;}
		unsigned int getIndex() const {return _idx;}

	private:
		unsigned int _idx;
		const bitset_pair* _status;
	};

	class const_pheno_iterator : public boost::iterator_facade<const_pheno_iterator, Phenotype const, boost::forward_traversal_tag>{

	public:
		const_pheno_iterator(unsigned int i, const PopulationManager& pop_mgr) : _idx(i), _mgr(pop_mgr), _itr(i, pop_mgr._pheno_status[i]){}

	private:
		friend class boost::iterator_core_access;

		void increment() { ++_idx;}
		bool equal(const const_pheno_iterator& other) const { return _idx == other._idx;}
		Phenotype const& dereference() const {_itr = Phenotype(_idx, _mgr._pheno_status[_idx]); return _itr;}

		unsigned int _idx;
		const PopulationManager& _mgr;
		mutable Phenotype _itr;

	};


	const_pheno_iterator beginPheno() const{
		return const_pheno_iterator(0, *this);
	}
	const_pheno_iterator endPheno() const{
		return const_pheno_iterator(_pheno_names.size(), *this);
	}


	explicit PopulationManager(const std::string& vcf_file);
	~PopulationManager(){}

	template <class T_cont>
	void loadLoci(T_cont& loci_out, const std::string& prefix, const std::string& sep, const Knowledge::Liftover::Converter* conv=0);

	// Usage functions
	unsigned int genotypeContribution(const Knowledge::Locus& locus) const;
	bool isRare(const Knowledge::Locus& locus, const bitset_pair& status, float lower, float upper) const;
	unsigned int getNumPhenotypes() const {return _pheno_names.size();}
	const std::string& getPhenotypeName(unsigned int i) const {return _pheno_names[i];}

	// Printing functions
	template <class Bin_cont>
	void printBins(std::ostream& os, const Bin_cont& bins, const Phenotype& pheno, const Knowledge::Information& info, const std::string& sep=",") const;
	template <class Bin_cont>
	void printBinsTranspose(std::ostream& os, const Bin_cont& bins, const Phenotype& pheno, const Knowledge::Information& info, const std::string& sep=",") const;

	static float c_phenotype_control;
	static std::string c_phenotype_file;

	// determine rarity of a variant by either case or control status
	static bool RareCaseControl;
	static bool NoSummary;

	static bool c_use_calc_weight;
	static bool _use_custom_weight;

	static float c_min_control_frac;

	static DiseaseModel c_model;
	static WeightModel c_weight_type;

private:

	void printEscapedString(std::ostream& os, const std::string& toPrint, const std::string& toRepl, const std::string& replStr) const;
	std::string getEscapeString(const std::string& sep) const;

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	// Loading functions
	unsigned int readVCFHeader(std::istream& v);
	void loadIndividuals();
	void parsePhenotypeFile(const std::string& filename);

	float getIndivContrib(const Knowledge::Locus& loc, int position, const bitset_pair& status, bool useWeights = false, const Knowledge::Information* const info = NULL, const Knowledge::Region* const reg = NULL) const;
	unsigned int getTotalContrib(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss=0) const;
	float getMAF(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss=0) const;
	float calcBrowningWeight(unsigned long N, unsigned long M) const;
	float calcWeight(const Knowledge::Locus& loc, const bitset_pair& status) const;
	float getCustomWeight(const Knowledge::Locus& loc, const Knowledge::Information& info, const Knowledge::Region* const reg = NULL) const;

	/*
	 * Fast atoi that handles up to 5 digits (max unsigned short is ~65K)
	 */
	unsigned short fast_atoi(const boost::string_ref& r){
		unsigned short value_ = 0;
		unsigned int len = r.size();
		switch(len){
        case  5:    value_ += (r[len- 5] - '0') * 10000;
        case  4:    value_ += (r[len- 4] - '0') * 1000;
        case  3:    value_ += (r[len- 3] - '0') * 100;
        case  2:    value_ += (r[len- 2] - '0') * 10;
        case  1:    value_ += (r[len- 1] - '0');
		}
		return value_;
	}

	boost::array<unsigned int, 2> getBinCapacity(Bin& bin, const bitset_pair& status) const;

	// mapping of ID -> position in the VCF file
	boost::unordered_map<std::string, int> _positions;
	// names of multiple phenotypes
	std::vector<std::string> _pheno_names;
	// the actual phenotypes as read from the phenotype file
	boost::unordered_map<std::string, std::vector<float> > _phenos;
	// phenotypes converted to case/control status
	std::vector<bitset_pair> _pheno_status;

	// the actual genotypes included in the VCF file
	// NOTE: this should eventually be a class that caches its results to a
	// temporary file for better use of memory!
	boost::unordered_map<const Knowledge::Locus*, bitset_pair > _genotypes;

	std::string _vcf_fn;
};

template <class Bin_cont>
void PopulationManager::printBinsTranspose(std::ostream& os, const Bin_cont& bins, const Phenotype& pheno, const Knowledge::Information& info, const std::string& sep) const{

	_use_custom_weight = (info.c_weight_files.size() > 0);
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

	boost::unordered_map<std::string, int>::const_iterator m_itr = _positions.begin();
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
	boost::unordered_map<std::string, std::vector<float> >::const_iterator pheno_status;
	while(m_itr != _positions.end()){
		float status = missing_status;
		pheno_status = _phenos.find((*m_itr).first);
		if (pheno_status != _phenos.end()){
			status = (*pheno_status).second[pheno.getIndex()];
		}
		os << sep << status;
		++m_itr;
	}

	os << "\n";

	typename Bin_cont::const_iterator b_itr = bins.begin();
	Bin::const_locus_iterator l_itr;
	boost::unordered_map<const Knowledge::Locus*, bitset_pair >::const_iterator loc_itr;

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

		// print for each person
		float bin_count = 0;
		m_itr = _positions.begin();
		while (m_itr != _positions.end()) {
			l_itr = (*b_itr)->variantBegin();
			bin_count = 0;
			while (l_itr != (*b_itr)->variantEnd()) {
				bin_count += getIndivContrib(**l_itr, (*m_itr).second, pheno.getStatus(), true, &info, (*b_itr)->getRegion());
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
void PopulationManager::printBins(std::ostream& os, const Bin_cont& bins, const Phenotype& pheno, const Knowledge::Information& info, const std::string& sep) const{
	static const float missing_status = std::numeric_limits<float>::quiet_NaN();
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

	boost::unordered_map<std::string, int>::const_iterator m_itr = _positions.begin();
	boost::unordered_map<std::string, int>::const_iterator m_end = _positions.end();

	boost::unordered_map<std::string, std::vector<float> >::const_iterator pheno_status;
	boost::unordered_map<std::string, std::vector<float> >::const_iterator pheno_end = _phenos.end();

	int pos;
	float status;
	float bin_count;

	while (m_itr != m_end){
		b_itr = bins.begin();
		b_end = bins.end();

		pos = (*m_itr).second;
		pheno_status = _phenos.find((*m_itr).first);
		status = missing_status;
		if (pheno_status != pheno_end){
			status = (*pheno_status).second[pheno.getIndex()];
		}

		printEscapedString(os, (*m_itr).first, sep, sep_repl);

		os << sep << status;

		while(b_itr != b_end){
			// Accumulate the contribution of this person in this bin
			l_itr = (*b_itr)->variantBegin();
			l_end = (*b_itr)->variantEnd();
			bin_count = 0;
			while(l_itr != l_end){
				bin_count += getIndivContrib(**l_itr, pos, pheno.getStatus(), true, &info, (*b_itr)->getRegion());
				++l_itr;
			}
			os << sep << std::setprecision(4) << bin_count;
			++b_itr;
		}

		os << "\n";
		++m_itr;
	}
}

template<class T_cont>
void PopulationManager::loadLoci(T_cont& loci_out, const std::string& prefix, const std::string& sep, const Knowledge::Liftover::Converter* conv){
	typedef std::string::const_iterator sc_iter;
	typedef boost::iterator_range<sc_iter> string_view;


	Utility::ICompressedFile vcf_f(_vcf_fn.c_str());
	unsigned int lineno = readVCFHeader(vcf_f);
	loadIndividuals();

	std::string geno_sep = "/";
	std::string alt_geno_sep = "|";

	std::string curr_line;
	unsigned int n_fields = 9 + _positions.size();

	std::string chr, id, ref, alt, filter, format;

	unsigned int bploc = 0;
	unsigned int gt_idx;
	unsigned int ft_idx;
	unsigned short curr_max;
	unsigned int max_count;
	bool lift_warn = false;
	std::ofstream unlift_out;
	std::vector<boost::string_ref> geno_list;
	std::vector<boost::string_ref> fields;
	std::vector<std::string> alleles;
	std::vector<boost::string_ref> call_list;
	std::vector<std::string> format_list;
	std::vector<std::pair<unsigned short, unsigned short> > calls;
	std::pair<unsigned short, unsigned short> curr_call;
	std::vector<unsigned int> call_count;


	const unsigned short missing_geno = static_cast<unsigned short>(-1);
	calls.reserve(_positions.size());
	format_list.reserve(10);
	call_list.reserve(2);

	//boost::char_separator<char> vcf_sep("\t");

	while(getline(vcf_f, curr_line)){
		++lineno;
		fields.clear();
		if(curr_line.size() > 2 && curr_line[0] != '#'){
			// In this case, we're looking at a marker

			boost::algorithm::iter_split(fields, curr_line, boost::first_finder("\t"));
			if(fields.size() != n_fields){
				std::cerr << "ERROR: Mismatched number of fields on line "
						<< lineno << std::endl;
				// throw exception here
				throw std::runtime_error("Mismatched number of fields in VCF file");
			}

			chr = std::string(fields[0].begin(), fields[0].end());
			bploc = boost::lexical_cast<unsigned int>(fields[1]);
			id = std::string(fields[2].begin(), fields[2].end());
			ref = std::string(fields[3].begin(), fields[3].end());
			alt = std::string(fields[4].begin(), fields[4].end());
			filter = std::string(fields[6].begin(), fields[6].end());
			format = std::string(fields[8].begin(), fields[8].end());

			// check for marker-level inclusion
			if(filter == "." || filter == "PASS"){

				// construct a locus object and lift over, if necessary
				Knowledge::Locus* loc = new Knowledge::Locus(chr,bploc,id, ref);
				if(conv){
					Knowledge::Locus* new_loc = conv->convertLocus(*loc);
					if (! new_loc){
						if(!lift_warn){
							std::string fn(prefix + "-unlifted.csv");
							unlift_out.open(fn.c_str());
							std::cerr << "WARNING: Some variants not lifted!  See "
									  << fn << " for details." << std::endl;

							unlift_out << "Chrom" << sep << "Pos" << sep
									   << "ID" << std::endl;
							lift_warn = true;

						}
						loc->print(unlift_out, sep);
						unlift_out << std::endl;

						delete loc;
						loc = 0;
					}else{
						delete loc;
						loc = new_loc;
					}
				}


				if(loc != 0){
					// get a list of all the alleles, in the correct order
					std::string allele_str = ref + "," + alt;
					boost::algorithm::split(alleles, allele_str, boost::is_any_of(","));
					// initialize the call count so I can easily determine the major allele
					call_count.clear();
					call_count.resize(alleles.size(), 0);

					// parse the format string
					format_list.clear();
					boost::algorithm::split(format_list, format, boost::is_any_of(":"));
					gt_idx = (find(format_list.begin(), format_list.end(), "GT") - format_list.begin());
					ft_idx = (find(format_list.begin(), format_list.end(), "FT") - format_list.begin());
					ft_idx = (ft_idx == format_list.size()) ? static_cast<unsigned int>(-1) : ft_idx;

					if(gt_idx == format_list.size()){
						std::cerr << "ERROR: No 'GT' format on line " << lineno <<
								", cannot continue." << std::endl;
						throw std::runtime_error("No GT given in format string");
					}

					calls.clear();
					for (unsigned int i=0; i<fields.size() - 9; i++){

						// parse the individual call
						geno_list.clear();
						//std::string currstr(fields[i+9].begin(), fields[i+9].end());
						boost::algorithm::iter_split(geno_list, fields[i+9], boost::first_finder(":"));

						if(fields[i+9] != "." && (ft_idx == static_cast<unsigned int>(-1) || geno_list[ft_idx] == "PASS")){
							call_list.clear();
							boost::algorithm::iter_split(call_list, geno_list[gt_idx], boost::first_finder(geno_sep));
							if(call_list.size() == 1){
								// we should be here very rarely!  If we're here,
								// we'll assume that the "primary" separator of
								// genotypes is in fact the "alternate", so swap them!
								boost::algorithm::iter_split(call_list, geno_list[gt_idx], boost::first_finder(alt_geno_sep));
								geno_sep.swap(alt_geno_sep);
							}

							if(call_list.size() != 2){
								std::cerr << "WARNING: Non-diploid found on line " <<
										lineno << ", setting to missing" << std::endl;
								calls.push_back(std::make_pair(missing_geno, missing_geno));
							} else {
								unsigned short g1, g2;
								if(call_list[0] != "." && call_list[1] != "."){
									g1 = fast_atoi(call_list[0]);
									g2 = fast_atoi(call_list[1]);
									calls.push_back(std::make_pair(g1, g2));
									++call_count[g1];
									++call_count[g2];
								} else{
									calls.push_back(std::make_pair(missing_geno, missing_geno));
								}
							}
						} else {
							calls.push_back(std::make_pair(missing_geno, missing_geno));
						}

					} // end iterating over genotypes

					// OK, now we'll find the major allele
					curr_max = 0;
					max_count = call_count[0];
					for(unsigned short i=1; i<call_count.size(); i++){
						if(call_count[i] > max_count){
							curr_max = i;
							max_count = call_count[i];
						}
					}

					// let's make sure that this isn't monoporphic
					std::sort(call_count.begin(), call_count.end());
					if(call_count.size() < 2 || call_count[call_count.size() - 2] == 0){
						delete loc;
					} else {

						bitset_pair curr_geno(std::make_pair(boost::dynamic_bitset<>(calls.size()), boost::dynamic_bitset<>(calls.size())));

						for(unsigned int i=0; i<calls.size(); i++){
							curr_call = calls[i];

							if(curr_call.first == missing_geno || curr_call.second == missing_geno){
								curr_geno.first.set(i);
								curr_geno.second.set(i);
							} else if (curr_call.first != curr_max && curr_call.second != curr_max){
								curr_geno.first.set(i);
							} else if (curr_call.first != curr_max || curr_call.second != curr_max) {
								curr_geno.second.set(i);
							}
						}

						// again, make sure it isn't monomorphic with regards to the
						// disease encoding
						if(getTotalContrib(curr_geno) == 0){
							delete loc;
						} else {
							loci_out.insert(loci_out.end(), loc);
							_genotypes.insert(std::make_pair(loc, curr_geno));
						}
					}
				}

			} // end marker parsing
		}

	} // end while(getline)

	if(lift_warn){
		unlift_out.close();
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
