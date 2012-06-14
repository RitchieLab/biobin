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

#include "binmanager.h"

#include "knowledge/Locus.h"

#include "vcftools/vcf_file.h"

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/dynamic_bitset.hpp>

using boost::dynamic_bitset;
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

	DataImporter(const string& filename) : vcf(filename, CompressedVCF) {}
	virtual ~DataImporter(){}
	
	typedef pair<dynamic_bitset<>, dynamic_bitset<> > bitset_pair;

	//bool open(const string& filename);
	
	/**
	 * Parses the VCF file for all variants.  Appends the variants to the input
	 * parameter given (allowing for multiple files to be read)
	 *
	 * @param loci_out The vector of Locus objects to append to
	 */
	template <class T_cont>
	void getLoci(T_cont& loci_out, const vector<bool>& controls=vector<bool>());

	static bool CompressedVCF;					///< gzipped file Y/N
	static bool KeepCommonLoci;
	// determine rarity of a variant by either case or control status
	static bool RareCaseControl;

	uint individualCount() {return vcf.N_indv;}
	//const vector<string>& getIndividualIDs() {return vcf.indv;}
	
	void parseSNP(Knowledge::Locus& loc, const dynamic_bitset<>& controls, bitset_pair& bitset_out, array<uint,2>& nonmissing_out);

	void remapLocus(Knowledge::Locus* orig_loc, Knowledge::Locus* new_loc);

private:

	// returns the minor allele frequency (frequency of the second most common
	// allele) given the counts and total count.
	float getMAF(const vector<int>& allele_count, uint nmcc);

	// No copying or assignment!
	DataImporter(const DataImporter& orig);
	DataImporter& operator=(const DataImporter& other);

	mutable VCF::vcf_file vcf;

	// A map to keep track of where in the file a locus resides
	unordered_map<Knowledge::Locus*, int> _locus_position;

};

template <class T_cont>
void DataImporter::getLoci(T_cont& loci_out, const vector<bool>& controls) {

	vector<bool> cases = controls;
	for(vector<bool>::iterator itr = cases.begin(); itr!= cases.end(); itr++){
		(*itr).flip();
	}

	set<string> unknownChromosomes;
	//T_cont::const_iterator pos = loci_out.end();

	uint totalSiteCount	= vcf.N_entries;

	//TODO: preallocate the map for some speed here
	//_locus_position.reserve(totalSiteCount);

	string line;
	vector<int> alleleCounts;
	uint nmcc = 0;					///< Just to avoid redundant conversions
	vector<int> alleleCounts_case;
	uint nmcc_case = 0;
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
			entry.get_allele_counts(alleleCounts_case, nmcc_case, cases, vcf.include_genotype[i]);
		}

		bool is_rare = (getMAF(alleleCounts, nmcc) < BinManager::mafCutoff) ||
				(RareCaseControl && nmcc_case > 0 && getMAF(alleleCounts_case, nmcc_case) < BinManager::mafCutoff);

		if(KeepCommonLoci || is_rare ){
			Locus* loc = new Locus(entry.get_CHROM(),entry.get_POS(),is_rare,entry.get_ID());

			loc->addAllele(entry.get_REF(), nmcc != 0 ? (alleleCounts[0] / ((double) nmcc) ) : 0);

			// If we are completely missing information (probably just for the controls),
			// Then get the allele information for the entire population, so
			// we can include the non-referent alleles with a frequency of -1
			if (nmcc == 0){
				entry.get_allele_counts(alleleCounts, nmcc, vcf.include_indv, vcf.include_genotype[i]);
				nmcc = 0;
			}

			//From here, they are all "ALT" alleles. We would have to evaluate an
			//if should we want to roll them all in together, since it's a different
			//call (with a different index scheme...)
			for (uint n = 1; n<alleleCount; n++)	{
				loc->addAllele(entry.get_ALT_allele(n-1),
						nmcc != 0 ? alleleCounts[n] / ((double) nmcc) : -1);
			}

			loci_out.insert(loci_out.end(), loc);
			_locus_position[loc] = i;
		}
	}

}

}


#endif	/* DATAIMPORTER_H */

