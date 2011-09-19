/* 
 * File:   binapplication.h
 * Author: torstees
 *
 * Created on June 22, 2011, 10:35 AM
 * 
 * One of the big TODO list things would be to integrate the two 
 * forms of SNPs: Biofilter and BioBin. Right now, we have two
 * very different approaches to SNPs. For the biofilter, we only
 * need a way to recognize names and associate them with a basepair
 * and chromosome. However, for biobin, we need to maintain alleles
 * and provide the ability to perform genotype conversion and 
 * some other stuff. So, the Locus object is much more complex. 
 * 
 * Most likely it's just a matter of moving the biobin Locus class
 * to someplace common and changing the biofilter code to use it...
 * 
 */

#ifndef BINAPPLICATION_H
#define	BINAPPLICATION_H

#include "binmanager.h"
#include "biofilter/application.h"
#include "dataimporter.h"
#include "individual.h"
#include <utility>

namespace BioBin {
	
class BinApplication : public Biofilter::Application {
public:
	BinApplication();
	virtual ~BinApplication();
	void SetReportPrefix(const char *pref);
	void InitVcfDataset(std::string& filename, 
			std::string& genomicBuild, 
			Knowledge::SnpDataset& lostSnps, 
			std::vector<uint>& locusRemap, 
			DataImporter& vcfimporter);
	//std::pair<uint, uint> _LoadVcfFile(std::string& filename, std::string& genomicBuild, Knowledge::SnpDataset& lostSnps);

	/**
	 * Return the individuals that have been loaded 
    * @return
    */
	const std::vector<Individual> &Individuals();

	/**
	 * Returns the number of SNPs that might contribute to a given bin
    * @param hits
    */
	void GetMaxBinHits(std::vector<uint>& hits);

	void GetBinContributors(std::vector<std::set<uint> >& contributors);

	/**
	 * Returns an array matching each of the bin names
    * @return
    */
	const Utility::StringArray& BinNames();

	/**
	 * Returns the region object at a given index
    * @param idx
    * @return
    */
	const Knowledge::Region& GetRegion(uint idx);

	Utility::StringArray phenotypeFilenames;			///< The file to be used to load phenotype values

	void ApplyPhenotypes();

	Utility::Locus &Locus(uint idx);

	/**
	 * Passes a list of genotypes (for all people at a given SNP) and returns the bin IDs 
    * @param snpIndex which SNP we are referring to
    * @param genotypes the original data from the vcf files
    * @param hits Vector containing the individual indexes where a variation occured
    * @param genotypes This is where we'll write genotype data
    * @return set of bin indexes for which this SNP applies
    */
	std::set<uint> ParseSNP(uint snpIndex, std::vector<char>& genotypes, std::vector<Individual>& data);
	
	/**
	 * Initialize the snp dataset based on the contents of the VCF file (this only sets up the locus portion, not individual data or bins)
    * @param filename
    * @param genomicBuild
    * @param lostSnps
    */
	std::pair<uint, uint> InitBins(std::vector<uint>& locusRemap, 
			DataImporter& vcfimporter);
	/**
	 * returns lookup region index -> snp index
	 */
	void GenerateBinContentLookup(std::multimap<uint, uint>& binContents);
	std::map<uint, uint> GetBinLookup();
private:
	BinManager binData;							///< Used to build and parse data into bins and genotypes
						///< Help to extract genotype data from vcf files
	std::vector<Individual> individuals;	///< This represents the actual data from the vcf files
	std::set<uint> binnable;					///< Just a map into the regions that created the entries
	std::map<uint, uint> binIndex;
	//std::vector<Utility::Locus> loci;		///< The loci associated with the dataset
	//Knowledge::SnpDataset loci;			/// This is now dataset
};

inline
void BinApplication::SetReportPrefix(const char *pref) {
	if (strcmp(pref, "") == 0)
		reportPrefix = "biobin";
	else
		reportPrefix = pref;
}

inline
void BinApplication::GetMaxBinHits(std::vector<uint>& hits) {
	std::vector<std::set<uint> > binContributors;
	binData.BuildContributorList(binContributors);
	uint binCount = binContributors.size();
	hits = std::vector<uint>(binCount, 0);
	
	for (uint i=0; i<binCount; i++) {
		hits[i] = binContributors[i].size();
	}
}

inline
Utility::Locus& BinApplication::Locus(uint idx) {
	return dataset[idx];
}
inline
std::map<uint, uint> BinApplication::GetBinLookup() {
	return binIndex;
}
inline
const std::vector<Individual>& BinApplication::Individuals() {
	return individuals;
}

inline
const Utility::StringArray& BinApplication::BinNames() {
	return binData.BinNames();
}

inline
const Knowledge::Region& BinApplication::GetRegion(uint idx) {
	return regions[idx];
}


inline
void BinApplication::GetBinContributors(std::vector<std::set<uint> >& contributors) {
	binData.BuildContributorList(contributors);
}

}
#endif	/* BINAPPLICATION_H */

