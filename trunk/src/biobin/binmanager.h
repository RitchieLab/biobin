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
#include <set>
#include <string>
#include <utility>

//#include "individual.h"
#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Group.h"
#include "knowledge/Locus.h"

#include "Bin.h"

//#include "knowledge/snpdataset.h"
//#include "knowledge/groupmanagerdb.h"
//#include "knowledge/regionmanagerdb.h"

using std::map;
using std::vector;
using std::string;
using std::pair;
using std::set;

namespace BioBin {

class Bin;

class BinManager {
public:
	typedef set<Bin*>::const_iterator const_iterator;

	BinManager();

	virtual ~BinManager();
	//BinManager(const BinManager& orig);
	
	void InitBins(const map<uint, Knowledge::GroupCollection*> &groups,
			const Knowledge::RegionCollection& regions,
			const vector<Knowledge::Locus*>& loci);

	int numRareVariants() const { return _rare_variants.size();}
	int numVariants() const {return _total_variants;}
	int numBins() const {return _bin_list.size();}
	//int numTotalVariants() const {return _total_variants;}

	const_iterator begin() const {return _bin_list.begin();}
	const_iterator end() const {return _bin_list.end();}

	void printBins(std::ostream& os, Knowledge::Locus* locus, const string& sep=":");



/*	void GenerateBins(std::map<Knowledge::Group*, Utility::IdCollection>& groups,
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
*/
	
	static uint IntergenicBinWidth;				///< The width of the intergenic bins within a chromosome
	static uint BinTraverseThreshold;			///< The number of SNPs to determine whether we continue traversing
	static uint MinBinSize;							///< How small do we tolerate bins to be before we ignore the bin altogether
	static bool ExpandByGenes;						///< Do we want to drop down to genes, if the group is large enough?
	static bool ExpandByExons;						///< Do we want to drop to to introns and exons, if the group is large enough
	static bool ExpandByFunction;					///< Indicate that we do want to use function 
	static float mafCutoff; 	///< Max maf to produce result in a bin

	// Carried this over from taskbincollapse, but I'm not sure what it does
	static uint maxSnpCount;
private:

	// Collapses all of the bins according to the preferences we set.
	void collapseBins();

	// The authoritative list of all of the bins.  Everything else holds pointers
	// to bins in this set.
	set<Bin*> _bin_list;
	// Mappint of Region IDs to bins
	map<int, Bin*> _region_bins;
	// Mapping of group IDs to bins
	map<int, Bin*> _group_bins;
	// List of intergenic bins
	map<pair<short, int>, Bin*> _intergenic_bins;
	// List of bins by locus
	map<Knowledge::Locus*, set<Bin*> > _locus_bins;

	// List of all rare variants
	set<Knowledge::Locus*> _rare_variants;
	int _total_variants;


/*	void CollectGroupLeaves(Knowledge::GroupManagerDB& gmgr,
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
*/
};

/*
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
*/
//inline

// * Sorts all of the variants found into either a "standard" variant or a "rare" variant.
// *
// * variants and rareVariants are sets of uints that are initially empty
/*
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
*/

} //namespace BioBin

#endif	/* BINMANAGER_H */

