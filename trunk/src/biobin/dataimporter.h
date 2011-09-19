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
#include "utility/strings.h"
#include "utility/types.h"
#include "vcftools/vcf_file.h"
#include <boost/algorithm/string.hpp>
#include "utility/locus.h"
namespace BioBin {
	


class DataImporter {
public:
	DataImporter();
	DataImporter(const DataImporter& orig);
	virtual ~DataImporter();
	
	/* Open/Close really just manages the datasource object
	 */
	bool Open(const char *filename, char chrom);
	bool Open(char chrom);
	void Close();
	
	/**
	 * Parses VCF file for allele frequencies
    * @param freqs
    * @return T/F for success
	 * 
	 * This is working under the assumption that we are doing all of this 
	 * for one chromosome at a time. If all chromosomes are in a single file
	 * then you can change chromosome and call these functions again...not
	 * most efficient, but it's how I am doing it! ;)
    */
	std::vector<Utility::Locus> GetAlleleFrequencies();
	void GetAllAlleleFrequencies(std::vector<Utility::Locus>& freq);
	/**
	 * returns array of count(Major alleles) for each locus
    * @param genotypes
    */
	void ParseSNP(uint snpIndex, std::vector<char>& genotypes);
	
	/**
	 * These are the biofilter style chromosomes...char indexes 1-25
	 */
	char chromosome;								///< Chromosome being parsed
	static bool CompressedVCF;					///< gzipped file Y/N
	void SetChromosome(char chrom);
	
	// Conversion for chromosomes['1','2',...'23','24','25'] 
	// If the user needs to define them in some other way...they should have a 
	// way to change this
	Utility::StringArray chromosomeNames;	
	
	/**
	 * Return the loci that have been scanned from this file
    * @return reference to the local locus vector
    */
	const std::vector<Utility::Locus>& GetLoci();
	
	
	uint IndividualCount();
	
	Utility::StringArray GetIndividualIDs();
private:
	uint totalIndividualEntries;				///< Number of individuals in the file(s)
	VCF::vcf_file *vcf;							///< This represents the vcf object we will be using
	VCF::vcf_entry *entry;						///< This is used to extract data out of the file

	void InitChromosomeNames();				///< Initialize the defaults
	std::vector<Utility::Locus> loci;		///< Help identify what is what within this region
	std::string filename;						///< Remember what file we are reading from
};

inline
DataImporter::DataImporter() : totalIndividualEntries(0), vcf(NULL), entry(NULL) { 
	InitChromosomeNames();
}

inline
DataImporter::DataImporter(const DataImporter& orig) : totalIndividualEntries(orig.totalIndividualEntries), vcf(orig.vcf), entry(orig.entry) {
	InitChromosomeNames();
}

inline
DataImporter::~DataImporter() { 
	if (entry)
		delete entry;
	if (vcf)
		delete vcf;
}

inline
Utility::StringArray DataImporter::GetIndividualIDs() {
	return vcf->indv;
}

inline
const std::vector<Utility::Locus>& DataImporter::GetLoci() {
	return loci;
}

inline
uint DataImporter::IndividualCount() {
	return totalIndividualEntries;
}

inline
void DataImporter::InitChromosomeNames() { 
	for (uint i=0;i<26; i++) {
		std::string s = Utility::ChromFromInt(i);
		boost::trim(s);
		chromosomeNames.push_back(s);
	}
}

inline
bool DataImporter::Open(char chromosome) {
	this->chromosome = chromosome;
	
	std::string chromosomeName = "";
	if (chromosome != -1)
		chromosomeName = chromosomeNames[chromosome];
	//I'm assuming we never want to exclude chromosomes...we'll just not parse them
	vcf = new VCF::vcf_file(std::string(filename), CompressedVCF, chromosomeName, "");
	totalIndividualEntries						= vcf->N_indv;
		
	entry = new VCF::vcf_entry(totalIndividualEntries);
	// VCF Tools allows the user to filter out individuals and SNPs based on certain 
	// criterion...we might want to do the same
	//vcf->apply_filters(params);
	
	return true;	
}

inline
bool DataImporter::Open(const char *filename, char chromosome) {
	if (entry)
		Close();
	this->filename = filename;
	return Open(chromosome);
}

inline
void DataImporter::Close() {
	if (entry) {
		delete entry;
		entry  = NULL;
	}
	
	if (vcf) {
		delete vcf;
		vcf = NULL;
	}
}

}


#endif	/* DATAIMPORTER_H */

