/* 
 * File:   binmanager.h
 * Author: torstees
 *
 * Created on July 21, 2011, 3:54 PM
 */

#ifndef BINMANAGER_H
#define	BINMANAGER_H

#include <vector>
#include <map>

#include "individual.h"
#include "knowledge/snpdataset.h"
#include "knowledge/groupmanagerdb.h"
#include "knowledge/regionmanagerdb.h"

namespace BioBin {

class BinManager {
public:
	BinManager();
	BinManager(const BinManager& orig);
	
/*	uint InitBins(std::map<uint, Knowledge::GroupManagerDB> &groups, 
			Knowledge::RegionManagerDB& regions, 
			Knowledge::SnpDataset& snps);
 */
	
	std::pair<uint, uint> InitBins(std::map<uint, Knowledge::GroupManagerDB> &groups, 
			Knowledge::RegionManagerDB& regions, 
			Knowledge::SnpDataset& snps,
			std::vector<uint> locusRemap);
	virtual ~BinManager();
	void GenerateBins(std::map<Knowledge::Group*, Utility::IdCollection>& groups, 
					std::map<uint, std::set<uint> > binnables, 
					std::vector<uint>& indexLookup, 
					uint &binID);
	void GenerateBins(std::map<uint, Utility::IdCollection>& binData, 
					Knowledge::RegionManagerDB& regions, 
					std::vector<uint>& indexLookup, 
					uint &binID, const std::string& type);

	void GenerateBin(Utility::IdCollection& binData, 
					std::vector<uint>& indexLookup, 
					uint &binID, 
					const std::string& name);
	void GenerateIntergenicBins(std::map<uint, std::map<uint, Utility::IdCollection> >& intergenic,  
				std::vector<uint>& indexLookup, uint& binID) ;
	std::set<uint> ParseSNP(uint snpIndex, std::vector<char>& genotypes, std::vector<Individual>& data);
	
	void BuildContributorList(std::vector<std::set<uint> >& contributors);
	
	void CollectVariantGroups(Utility::IdCollection& variants, Utility::IdCollection& rareVariants);
	
	void DescribeLocus(uint snpIndex, std::ostream& os, Knowledge::RegionManagerDB& regions, Knowledge::SnpDataset& snps);
	
	//std::set<uint> ParseSNP(uint snpIndex, std::vector<char>& genotypes, std::vector<Individual>& data);
	
	const Utility::StringArray& BinNames();
	std::string BinName(uint index);

	const std::vector<uint>& ImplicationIndex();
	
	
	static uint IntergenicBinWidth;				///< The width of the intergenic bins within a chromosome
	static uint BinTraverseThreshold;			///< The number of SNPs to determine whether we continue traversing
	static uint MinBinSize;							///< How small do we tolerate bins to be before we ignore the bin altogether
	static bool ExpandByGenes;						///< Do we want to drop down to genes, if the group is large enough?
	static bool ExpandByExons;						///< Do we want to drop to to introns and exons, if the group is large enough
	static bool ExpandByFunction;					///< Indicate that we do want to use function 
	static float mafCutoff;							///< Max maf to produce result in a bin
private:
	void CollectGroupLeaves(Knowledge::GroupManagerDB& gmgr,
		std::map<uint, Utility::IdCollection>& regionLookup,
		std::map<uint, std::set<uint> >& regionToBinnable,
		std::map<Knowledge::Group*, Utility::IdCollection>& leaves,
		Utility::IdCollection& genesUsed);

	void CollectGroupLeaves(Knowledge::GroupManagerDB& gmgr,
		std::map<uint, Utility::IdCollection>& regionLookup,
		std::map<uint, std::set<uint> >& regionToBinnable,
		std::map<Knowledge::Group*, Utility::IdCollection>& leaves,
		Utility::IdCollection& genesUsed,
		Utility::IdCollection& visited, 
		uint groupIdx);

	void BinSNPs(Utility::IdCollection& snpIdx, 
		Knowledge::SnpDataset& snps, 
		std::map<uint, uint> locusRemap, 
		std::string& name,
		uint &binID);
	/**
	 * All but a small few will point to only one bin
	 */
	std::map<uint, uint> genotypeMap;			///< snp index to genotype index
	std::vector<uint> binIDs;						///< The general bin
	std::vector<uint> sourceIndex;				///< where to get the raw data
	std::multimap<uint, uint> multiBins;		///< Handle rare SNPs in more than one bin
	Utility::StringArray binNames;
	std::vector<uint> implicationIndex;			///< Bin's implication score
};


inline
void BinManager::BuildContributorList(std::vector<std::set<uint> >& contributors) {
	uint binCount = binNames.size();
	contributors.clear();
	contributors.resize(binCount);
	std::map<uint, uint>::iterator notMulti = multiBins.end();
	
	uint genotype = (uint)-1;
	uint snpCount = binIDs.size();
	for (uint i=0; i<snpCount; i++) {
		uint &id = binIDs[i];
		if (id != genotype) {
			if (multiBins.find(i) == notMulti)
				contributors[id].insert(i);
		}
	}
	
	std::map<uint, uint>::iterator itr = multiBins.begin();
	while (itr != notMulti) {
		contributors[itr->second].insert(itr->first);
		itr++;
	}
}

inline
const std::vector<uint>& BinManager::ImplicationIndex() {
	return implicationIndex;
}

inline
std::string BinManager::BinName(uint index) {
	return binNames[index];
}

inline
const Utility::StringArray& BinManager::BinNames() {
	return binNames;
}

inline
/**
 * Sorts all of the variants found into either a "standard" variant or a "rare" variant.
 *
 * variants and rareVariants are sets of uints that are initially empty
 */
void BinManager::CollectVariantGroups(Utility::IdCollection& variants, Utility::IdCollection& rareVariants) {
	uint count = binIDs.size(); 
	uint notRare = (uint)-1;
	
	for (uint i=0; i<count; i++) {
		std::cerr<<"asdfasdf "<<binIDs[i]<<"\t"<<variants.size()<<"\t"<<rareVariants.size()<<"\n";
		if (binIDs[i] == notRare)
			variants.insert(i);
		else
			rareVariants.insert(i);
	}
	
}

}

#endif	/* BINMANAGER_H */

