/* 
 * File:   SnpDataset.h
 * Author: torstees
 *
 * Provides compact SNP storage and flexible snp queries. All external SNP references
 * involve SNP index (currently an unsigned integer...so ~4 billion) When
 * something needs the actual SNP record, it will query for it using one of the
 * SNP queries.
 * 
 * Created on March 3, 2011, 9:17 AM
 */

#ifndef SNPDATASET_H
#define	SNPDATASET_H

#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include "utility/locus.h"
//#include "utility/strings.h"
#include <sstream>
#include "liftover/snp.h"
#include "chromosome.h"

#define VARIATION_FORMAT_VERSION 1.0


namespace Knowledge {


class SnpDataset {
public:
	typedef std::map<char, Utility::IdCollection> ChromBasedIdCollection;
	static bool BinaryArchive;

	SnpDataset(const char *variationsFN = "");
	SnpDataset(const SnpDataset& other);

	/**
	 * LoadData builds the SNP array. This is required prior to any lookups
    * @param rsids
	 * @param snpsFound These are the RS Numbers (literally, numbers) that were found in the variations file.
	 * @note There is no way to distinguish between RS numbers which appear multiple times....so they all will be stored
    * @return
    */
	uint LoadData(const std::set<std::string>& rsids, std::set<std::string>& snpsFound);
	uint LoadData(const char *filename, const std::set<std::string>& rsids, std::set<std::string>& snpsFound);
	uint LoadMapData(const char *filename, uint uscBuildVersion, bool performAlignment);
	/**
	 * LoadData from a SnpArray-which is likely the result of parsing a marker info file
	 * The system doesn't actually need to do anything except build out the internal
	 * representation and uplift the positions to match the build in the database
    * @param data
    * @param uscBuildVersion
    * @return
    */
	uint LoadData(const Utility::SnpArray& data, uint uscBuildVersion);

	void LoadFromString(const char *input);

	/**
	 * Store the SNP data within a text/binary map file
    * @param filename
    * @param binary T/F the file is stored as binary
    */
	void LoadArchive(const char *filename);


	/**
	 * Load SNP data from text/binary map file
    * @param filename
	 * @param uses english representation of the role if true
    * @param binary T/F the file is binary
    */
	void SaveArchive(const char *filename, bool writeRole = false);

	/** 
	 * (Assumes pared variations file) returns set of indexes
    * @param rsids
    * @param snps
    * @return how many were found
    */
	void RsToSnpIndexes(const std::set<std::string>& rsids, Utility::IdCollection& snps);

	/**
	 * Return a list of indexes that fall between begin and end on a given chromosome
    * @param chrom
    * @param begin
    * @param end
    * @param snps	This is not cleared, so it can be used in an additive fashion
    */
	void RangeSnpLookup(char chrom, uint begin, uint end, Utility::IdCollection& snps) const;
	void RangeSnpLookup(char chrom, uint begin, uint end, std::set<LiftOver::SNP>& snps) const;
	void PositionLookup(Utility::IdCollection& snps, Utility::IdCollection& pos) const;
	void RoleDescription(uint id, const char *desc);
	std::string RoleDescription(uint id);

	/**
	 * Easy way of realizing a SNP based on it's index
    * @param index
    * @return
	 * @note For now, I'm not checking the bounds...it's up to the client to know what he/she is doing. This is not java, it's efficient:P
    */
	Utility::Locus& operator[](uint index);

	/**
	 * Empties all elements from the dataset
    */
	void Clear();

	/**
	 * Returns the number of SNPs in the dataset
    * @return
    */
	uint Size() const;

	virtual ~SnpDataset() {}

	uint AddSNP(const Utility::Locus& l);
	uint AddSNP(char chrom, uint pos, const char *rsid, int role = 0);
	void AddSNPs(const Utility::SnpArray& snps);
	
	/**
	 * Reconciles the uplifted data when possibly by scanning the SNP data from the variations file.
	 * SNPs are categorized into 4 groups:
	 *		Success		- These are not reported
	 *		Off by one	- These have only a single bp difference from the variance information
	 *						- These SNPs are updated to the variance locations and counted under OBO
	 *		Tolerant		- These are off by more than 1, but within tolerance.
	 *						- These will be counted under tolerant
	 *		Failures		- These are SNPs whose location differs by more than tolerance
	 *						- These SNPs are removed from the dataset and always reported by Chrom/RS of original location
	 *		Unknown		- These are SNPs which are not found in the variation file. 
	 *						- These remain unchanged from the liftover results
	 *		
	 */
	uint ReconcileLiftover(std::multimap<Utility::Locus, Utility::Locus>& converted, std::ostream& os, int tolerance);
	
	/**
	 * Attempts to identify the role information from SNPs added via AddSNP
    * @param errorWidth is the amount of difference allowed between the same RS Number and the one listed in the map file
    * @return 
    */
	uint AlignData();

	static uint rsMapToPositionTolerance;
	static bool detailedReport;
	void WriteMarkerInfo(const char *filename, char sep = '\t');
	void SetVariationsFilename(const char *variationsFN);

	void GetFragment(const std::set<std::string>& rsids, SnpDataset& other);
	
	/**
	 * Associate a region with on the chromosome with an index
    * @param begin
    * @param end
    * @param regionIdx
    */
	void AddRegion(char chrom, uint begin, uint end, uint regionIdx);
	
	bool GetRegionCoverage(uint chrom, uint point, Utility::IdCollection& regionIdxs);
	
	bool GetRegionCoverage(uint snpIndex, Utility::IdCollection& regionIdxs);
	void InitTerms();
protected:
	/**
	 * Used to convert the integer from the variations file into an RSID
    * @param id - this is assumed to be the RS ID without the letters RS
    * @return 
    */
	std::string MakeRSID(uint id);
	uint ExtractRsInt(const std::string& rsid);
	uint FindSnp(char chrom, const char *rsid, uint pos);
	uint FindSnp(char chrom, const char *rsid, uint pos, uint tolerance);
	/**
	 * Store the SNP data within a text/binary map file
    * @param filename
    * @param binary T/F the file is stored as binary
    */
	void LoadArchiveBinary(const char *filename);

	/**
	 * Load SNP data from text/binary map file
    * @param filename
    * @param binary T/F the file is binary
    */
	void SaveArchiveBinary(const char *filename);

	std::string variationsFilename;
	Chromosome chromosomes[Utility::Locus::MaxChromosomes];
	std::map<uint, std::string> roleDescription;
	Utility::SnpArray markers;
	uint fileVersion;
	std::set<std::string> exonTerms;
	std::set<std::string> intronTerms;
	std::set<std::string> regulatoryTerms;
};

inline
SnpDataset::SnpDataset(const SnpDataset& other)
		: variationsFilename(other.variationsFilename),
		  chromosomes(other.chromosomes),
		  roleDescription(other.roleDescription),
		  markers(other.markers),
		  fileVersion(other.fileVersion) {}

/****************************** SNP **********************************/


inline
void SnpDataset::InitTerms() {
	std::string eTerms			= "STOP_GAINED,STOP_LOST,COMPLEX_INDEL,NON_SYNONYMOUS_CODING,FRAMESHIFT_CODING,PARTIAL_CODON,SYNONYMOUS_CODING,CODING_UNKNOWN,5PRIME_UTR,3PRIME_UTR";
	std::string iTerms			= "ESSENTIAL_SPLICE_SITE,INTRONIC,WITHIN_NON_CODING_GENE";
	std::string rTerms			= "WITHIN_MATURE_miRNA,NMD_TRANSCRIPT,UPSTREAM,DOWNSTREAM";
	
	exonTerms						= Utility::ToSet<std::string>(eTerms.c_str(), ",");
	intronTerms						= Utility::ToSet<std::string>(iTerms.c_str(), ",");
	regulatoryTerms				= Utility::ToSet<std::string>(rTerms.c_str(), ",");
}



inline
std::string SnpDataset::RoleDescription(uint id) {
	return roleDescription[id];
}

inline
std::string SnpDataset::MakeRSID(uint id) {
	std::stringstream ss;
	ss<<"RS"<<id;
	return ss.str();
}

inline
bool SnpDataset::GetRegionCoverage(uint snpIndex, Utility::IdCollection& regionIdxs) {
	return GetRegionCoverage(markers[snpIndex].chrom, markers[snpIndex].pos, regionIdxs);
}

inline
uint SnpDataset::ExtractRsInt(const std::string& rsid) {
	std::string rs = rsid;
	if (rs.find("r") != std::string::npos || rs.find("R") != std::string::npos)
			rs.erase(0,2);
	uint id = atoi(rs.c_str());
	assert(id);
	return id;

}
}
#endif	/* SNPDATASET_H */

