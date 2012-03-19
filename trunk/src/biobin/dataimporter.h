/* 
 * File:   dataimporter.h
 * Author: torstees
 *
 * This is a simple interface into the vcf tools that will allow 
 * the application to load allele frequencies and allelic information for
 * vcf files. 
 * 
 * Each importer will have a single file, and is probably going to be tied
 * to a certain chromosome. There really shouldn't be anything stored other
 * than a vcf_file object. 
 * 
 * Caveats:
 * There is an issue with how I'm storing data that could make for trouble.
 * Genotypes and bins are built and populated based on the local data and 
 * not upon the referent allele. So, if the user is planning on running
 * analysis on data from two separate VCF files, they will likely get
 * some bad output where the two files different in allele frequency: (i.e.
 * bin layout will probably change between the two and the way genotypes are
 * called could also change, depending on which is determined to be the 
 * major allele (not to mention the fact that allele frequency determines 
 * which SNPs are considered to be genotypes)
 * 
 * 
 * Created on June 16, 2011, 4:21 PM
 */

#ifndef DATAIMPORTER_H
#define	DATAIMPORTER_H

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "knowledge/Locus.h"

#include "vcftools/vcf_file.h"

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>

using std::string;
using std::vector;
using std::map;
using std::pair;
using boost::unordered_map;
using boost::array;

using Knowledge::Locus;

namespace BioBin {
	
class DataImporter {
public:

	DataImporter(const string& filename) : vcf(filename) {}
	virtual ~DataImporter(){}
	
	//bool open(const string& filename);
	
	/**
	 * Parses the VCF file for all variants.  Appends the variants to the input
	 * parameter given (allowing for multiple files to be read)
	 *
	 * @param loci_out The vector of Locus objects to append to
	 */
	template <class T_cont>
	void getLoci(T_cont& loci_out, const vector<bool>& controls=vector<bool>());

	template <class T_cont>
	void getCaseAF(const T_cont& loci, const vector<bool>& controls,
			unordered_map<Knowledge::Locus*, float>& maf_out) const;

	template <class T_cont>
	void getNumNonMissing(const T_cont& loci, const vector<bool>&controls,
			unordered_map<Knowledge::Locus*, array<uint,2> >& num_out);

	static bool CompressedVCF;					///< gzipped file Y/N
	static bool KeepCommonLoci;

	uint individualCount() {return vcf.N_indv;}
	const vector<string>& getIndividualIDs() {return vcf.indv;}
	
	void parseSNP(Knowledge::Locus& loc, vector<short>& genotypes_out);

private:

	// No copying or assignment!
	DataImporter(const DataImporter& orig);
	DataImporter& operator=(const DataImporter& other);

	mutable VCF::vcf_file vcf;

	// A map to keep track of where in the file a locus resides
	unordered_map<Knowledge::Locus*, int> _locus_position;

};

template <class T_cont>
void DataImporter::getLoci(T_cont& loci_out, const vector<bool>& controls) {



	set<string> unknownChromosomes;
	//T_cont::const_iterator pos = loci_out.end();

	uint totalSiteCount	= vcf.N_entries;

	//TODO: preallocate the map for some speed here
	//_locus_position.reserve(totalSiteCount);

	string line;
	vector<int> alleleCounts;
	double nonMissingChrCount = 0.0;
	uint nmcc = 0;					///< Just to avoid redundant conversions
	VCF::vcf_entry entry(vcf.N_indv);

	for (uint i=0; i<totalSiteCount; i++) {
		vcf.get_vcf_entry(i, line);
		entry.reset(line);
		entry.parse_basic_entry(true);
		entry.parse_genotype_entries(true);
		uint alleleCount = entry.get_N_alleles();

		// To deal with the default parameter, we really just want to use
		// vcf.include_indivs, but NOOOO C++ has to be a pain
		if(controls.size() == 0){
			entry.get_allele_counts(alleleCounts, nmcc, vcf.include_indv, vcf.include_genotype[i]);
		}else{
			entry.get_allele_counts(alleleCounts, nmcc, controls, vcf.include_genotype[i]);
		}
		nonMissingChrCount = nmcc;

		Locus* loc = new Locus(entry.get_CHROM(),entry.get_POS(),entry.get_ID());

		//if (chr == 0) // TODO Determine how to handle these that we don't recognize. We need to avoid pulling them when we pull genotypes
		//	unknownChromosomes.insert(entry->get_CHROM());

		loc->addAllele(entry.get_REF(), nmcc != 0 ? (alleleCounts[0] / nonMissingChrCount) : 0);
		if (nmcc == 0){
			entry.get_allele_counts(alleleCounts, nmcc, vcf.include_indv, vcf.include_genotype[i]);
			nmcc = 0;
		}

		//From here, they are all "ALT" alleles. We would have to evaluate an
		//if should we want to roll them all in together, since it's a different
		//call (with a different index scheme...)
		for (uint n = 1; n<alleleCount; n++)	{
			loc->addAllele(
					entry.get_ALT_allele(n-1),
					nmcc != 0 ? alleleCounts[n] / nonMissingChrCount : -1);
		}

		loci_out.insert(loci_out.end(), loc);
		_locus_position[loc] = i;
	}

}

template <class T_cont>
void DataImporter::getCaseAF(const T_cont& loci, const vector<bool>& controls,
		unordered_map<Knowledge::Locus*, float>& maf_out) const{
	int num_cases = 0;
	vector<bool> cases = controls;

	for (unsigned int i=0; i < controls.size(); i++){
		cases[i].flip();
		num_cases += cases[i];
	}
	if(num_cases == 0){
		std::cerr << "WARNING: No cases found!  "
				"No data to calculate case allele frequency" << std::endl;
	}

	typename T_cont::const_iterator l_itr = loci.begin();
	typename T_cont::const_iterator l_end = loci.end();

	vector<int> alleleCounts;
	string line;
	double nonMissingChrCount = 0.0;
	uint nmcc = 0;					///< Just to avoid redundant conversions
	VCF::vcf_entry entry(vcf.N_indv);
	float missing_val = -1;

	unordered_map<Knowledge::Locus*, int>::const_iterator locus_pos_itr;
	unordered_map<Knowledge::Locus*, int>::const_iterator locus_pos_end =
			_locus_position.end();

	while(l_itr != l_end){

		locus_pos_itr = _locus_position.find(*l_itr);
		if(locus_pos_itr == locus_pos_end){
			std::cerr << "WARNING: Could not find " << (*l_itr)->getID() <<
					" when calculating case AF" << std::endl;

			maf_out[*l_itr] = missing_val;
		} else {

			vcf.get_vcf_entry((*locus_pos_itr).second, line);
			entry.reset(line);
			entry.parse_basic_entry(true);
			entry.parse_genotype_entries(true);

			entry.get_allele_counts(alleleCounts, nmcc, cases,
					vcf.include_genotype[(*locus_pos_itr).second]);
			nonMissingChrCount = nmcc;

			if(nmcc > 0){
				maf_out[*l_itr] =  1 - (alleleCounts[(*l_itr)->getMajorPos()] / nonMissingChrCount);
			}else{
				maf_out[*l_itr] = missing_val;
			}
		}
		++l_itr;
	}
}

template <class T_cont>
void DataImporter::getNumNonMissing(const T_cont& loci, const vector<bool>&controls,
		unordered_map<Knowledge::Locus*, array<uint, 2> >& num_out){

	vector<bool> cases = controls;

	for (unsigned int i=0; i < controls.size(); i++){
		cases[i].flip();
	}

	typename T_cont::const_iterator l_itr = loci.begin();
	typename T_cont::const_iterator l_end = loci.end();

	vector<int> alleleCounts;
	string line;
	uint nmcc_case = 0;					///< Just to avoid redundant conversions
	uint nmcc_cont  = 0;
	VCF::vcf_entry entry(vcf.N_indv);

	unordered_map<Knowledge::Locus*, int>::const_iterator locus_pos_itr;
	unordered_map<Knowledge::Locus*, int>::const_iterator locus_pos_end =
			_locus_position.end();

	array<uint,2> case_cont_array;

	int num_not_found = 0;

	while(l_itr != l_end){

		locus_pos_itr = _locus_position.find(*l_itr);
		if(locus_pos_itr == locus_pos_end){
			std::cerr << "WARNING: Could not find " << (*l_itr)->getID() <<
					" when calculating case AF" << std::endl;

			//should create array of zeros
			num_out[*l_itr];
		} else {

			vcf.get_vcf_entry((*locus_pos_itr).second, line);
			entry.reset(line);
			entry.parse_basic_entry(true);
			entry.parse_genotype_entries(true);

			entry.get_allele_counts(alleleCounts, nmcc_case, cases,
					vcf.include_genotype[(*locus_pos_itr).second]);
			entry.get_allele_counts(alleleCounts, nmcc_cont, controls,
					vcf.include_genotype[(*locus_pos_itr).second]);

			case_cont_array[0] = nmcc_cont;
			case_cont_array[1] = nmcc_case;
			if(case_cont_array[0] == 0 || case_cont_array[1] == 0){
				++num_not_found;
			}
			num_out[*l_itr] = case_cont_array;
		}
		++l_itr;
	}
	if(num_not_found){
		std::cerr << "WARNING: Data completely missing for cases or controls "
				<< "for " << num_not_found << " loci.  See the locus.csv report"
				<< " for details.\n";
	}
}

}


#endif	/* DATAIMPORTER_H */

