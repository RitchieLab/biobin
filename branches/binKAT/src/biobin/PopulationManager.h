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
#include <boost/thread/mutex.hpp>

#include "Bin.h"

#include "util/ICompressedFile.h"
#include "util/string_ref.hpp"
#include "util/Phenotype.h"

#include "knowledge/Locus.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/Information.h"
#include "knowledge/Region.h"

namespace BioBin{

class BinManager;

namespace Test{
class Test;
}

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

	class const_pheno_iterator: public boost::iterator_facade<
			const_pheno_iterator, Utility::Phenotype const,
			boost::forward_traversal_tag> {

	public:
		const_pheno_iterator(unsigned int i, const PopulationManager& pop_mgr) :
			_idx(i), _mgr(pop_mgr), _itr(i, pop_mgr._pheno_status[i]) {
		}

	private:
		friend class boost::iterator_core_access;

		void increment() {
			++_idx;
		}
		bool equal(const const_pheno_iterator& other) const {
			return _idx == other._idx;
		}
		Utility::Phenotype const& dereference() const {
			_itr = Utility::Phenotype(_idx, _mgr._pheno_status[_idx]);
			return _itr;
		}

		unsigned int _idx;
		const PopulationManager& _mgr;
		mutable Utility::Phenotype _itr;

	};


	typedef std::vector<std::string>::const_iterator const_sample_iterator;
	const_sample_iterator beginSample() const{
		return _sample_names.begin();
	}

	const_sample_iterator endSample() const{
		return _sample_names.end();
	}

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
	unsigned short getIndivGeno(const Knowledge::Locus& loc, int position) const;
	float getAvgGenotype(const Knowledge::Locus& locus, const boost::dynamic_bitset<>* nonmiss_status=0) const;
	bool isRare(const Knowledge::Locus& locus, const bitset_pair& status, float lower, float upper) const;
	unsigned int getNumPhenotypes() const {return _pheno_names.size();}
	unsigned int getNumCovars() const {return _covar_names.size();}
	unsigned int getNumSamples() const {return _sample_names.size();}
	unsigned int getSamplePosition(const std::string& s) const;
	const std::string& getPhenotypeName(unsigned int i) const {return _pheno_names[i];}
	float getPhenotypeVal(const std::string& sample, const Utility::Phenotype& pheno) const;
	const std::vector<float>& getCovariates(const std::string& sample) const;
	float getTotalIndivContrib(const Bin& b, int pos, const Utility::Phenotype& pheno) const;
	float getLocusWeight(const Knowledge::Locus& loc, const Utility::Phenotype& pheno, const Knowledge::Region* reg=NULL) const;

	// working with the Knowlede::Information
	void setInfo(const Knowledge::Information* info) {
		_info = info;
		_use_custom_weight = _info && _info->c_weight_files.size() > 0;
	}
	const Knowledge::Information* getInfo() const { return _info;}

	// Printing functions
	void printBins(std::ostream& os, const BinManager& bins, const Utility::Phenotype& pheno, const std::string& sep=",") const;
	void printBinsTranspose(std::ostream& os, const BinManager& bins, const Utility::Phenotype& pheno, const std::string& sep=",") const;

	static float c_phenotype_control;
	static std::string c_phenotype_file;
	static std::string c_covariate_file;

	// determine rarity of a variant by either case or control status
	static bool RareCaseControl;
	static bool NoSummary;

	static bool c_keep_monomorphic;
	static bool c_use_calc_weight;

	static float c_min_control_frac;

	static DiseaseModel c_model;
	static WeightModel c_weight_type;

	static std::vector<const Test::Test*> c_tests;

	const static unsigned short missing_geno = static_cast<unsigned short>(-1);

private:

	void printEscapedString(std::ostream& os, const std::string& toPrint, const std::string& toRepl, const std::string& replStr) const;
	std::string getEscapeString(const std::string& sep) const;

	// NO copying or assignment!
	PopulationManager(const PopulationManager&);
	PopulationManager& operator=(const PopulationManager&);

	// Loading functions
	unsigned int readVCFHeader(std::istream& v);
	void loadIndividuals();
	void parseTraitFile(const std::string& fn,
			std::vector<std::string>& names_out,
			boost::unordered_map<std::string, std::vector<float> >& vals_out,
			const std::string& var_prefix="pheno");

	float getIndivContrib(const Knowledge::Locus& loc, int position, const Utility::Phenotype& pheno, bool useWeights = false, const Knowledge::Region* const reg = NULL) const;
	unsigned int getTotalContrib(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss=0) const;
	float getMAF(const bitset_pair& geno, const boost::dynamic_bitset<>* nonmiss=0) const;
	float calcBrowningWeight(unsigned long N, unsigned long M) const;
	float calcWeight(const Knowledge::Locus& loc, const bitset_pair& status) const;
	float getCustomWeight(const Knowledge::Locus& loc, const Knowledge::Region* const reg = NULL) const;

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

	// Note: thefollowing 2 variables are inverses of each other, so:
	// i == _positions[_sample_names[i]]
	// s == _sample_names[_positions[s]]

	// mapping of ID -> position in the VCF file
	boost::unordered_map<std::string, unsigned int> _positions;
	// vector of samples, as given in the VCF file
	std::vector<std::string> _sample_names;

	// names of multiple phenotypes
	std::vector<std::string> _pheno_names;
	// the actual phenotypes as read from the phenotype file
	boost::unordered_map<std::string, std::vector<float> > _phenos;
	// phenotypes converted to case/control status
	std::vector<bitset_pair> _pheno_status;

	// names of the covariates
	std::vector<std::string> _covar_names;
	//unsigned int _num_covars;
	// covariates as read in the covariate file(s)
	boost::unordered_map<std::string, std::vector<float> > _covars;

	// the actual genotypes included in the VCF file
	// NOTE: this should eventually be a class that caches its results to a
	// temporary file for better use of memory!
	boost::unordered_map<const Knowledge::Locus*, bitset_pair > _genotypes;

	std::string _vcf_fn;

	bool _use_custom_weight;

	const Knowledge::Information* _info;

};

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

	while(vcf_f.good()){
		getline(vcf_f, curr_line);
		++lineno;
		fields.clear();
		if(curr_line.size() > 2 && curr_line[0] != '#'){
			// In this case, we're looking at a marker

			boost::algorithm::iter_split(fields, curr_line, boost::first_finder("\t"));
			if(fields.size() != n_fields){
				std::cerr << "ERROR: Mismatched number of fields on line "
						<< lineno << std::endl;
				std::cerr << "Expected # of fields: "<<n_fields << std::endl;
				std::cerr << "Seen # of fields: " << fields.size() << std::endl;
				std::cerr << "Line size: " << curr_line.size() << std::endl;
				std::cerr << "Last field: " << fields[fields.size() - 1] << std::endl;
				std::cerr << "2nd to Last field: " << fields[fields.size() - 2] << std::endl;
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
				Knowledge::Locus* loc = new Knowledge::Locus(chr,bploc,id,ref);
				if(conv){
					Knowledge::Locus* new_loc = conv->convertLocus(*loc);
					// make sure to drop loci that lift to unknown chromosomes, too!
					if (! new_loc || new_loc->getChrom() == Knowledge::LocusPosition::UNKNOWN_CHR_RETURN){
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
					if(!c_keep_monomorphic && (call_count.size() < 2 || call_count[call_count.size() - 2] == 0)){
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
						if(!c_keep_monomorphic && getTotalContrib(curr_geno) == 0){
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
