/* 
 * File:   binmanager.cpp
 * Author: torstees
 * 
 * Created on July 21, 2011, 3:54 PM
 */

#include "binmanager.h"

namespace BioBin {

uint BinManager::IntergenicBinWidth = 50000;
uint BinManager::BinTraverseThreshold = 50;
uint BinManager::MinBinSize = 1;
bool BinManager::ExpandByExons = true;
bool BinManager::ExpandByFunction = true;
float BinManager::mafCutoff = 0.05;

BinManager::BinManager() {
}

BinManager::BinManager(const BinManager& orig) {
}

BinManager::~BinManager() {
}

/**
 * Basically, we'll need to generate a list of "leaf" nodes for
 * the selected set of groups. We'll store these as pointers to 
 * the group in a set or a vector.
 * 
 * If the user wishes to drill down further (to the gene level), 
 * we'll walk through the list of "leaves" and check each one
 * for how many snps fall inside each. If it is larger than
 * the threshold, we'll replace it with it's geneBins for each
 * of the genes. 
 * 
 * If the user wishes to drill further yet, we'll check each
 * of the genes for the presence of 2 or more subgroups and
 * replace them with their intronic/exonic oriented bins.
 * 
 * Finally, if the user wishes to express by functionality, 
 * each bin will be evaluated for the number of function based
 * bins and the bin entry will be replaced by 2 or more
 * functional bins
 * 
 * At any time of expansion/expression, we'll just drop a 
 * bin if it's contents falls below the suggested minimum
 * threshold.
 * 
 * Once the bin list(s) are created, we'll simply construct
 * the snpIndex->bin lookup and destroy the intermediate bin
 * stuff. 
 * 
 * Locus Remap shows what the current index is and what it was. 
 * That is the index we actually store the index->bin/genotype
 * information for
 */
std::pair<uint, uint> BinManager::InitBins(
		std::map<uint, Knowledge::GroupManagerDB> &groups,
		Knowledge::RegionManagerDB& regions,
		Knowledge::SnpDataset& snps,
		std::vector<uint> locusRemap) {
	/**
	 * The following items will result in a bin once we reach the end of the 
	 * function. So, if we refine an item at one level, we should remove it from
	 * the original structure. These will be the index into the dataset-not into
	 * the vcf file (the locusRemap is used at the very last step only)
	 */
	/** Group -> snpIdx*/
	std::map<Knowledge::Group*, Utility::IdCollection> leaves;
	std::map<uint, Utility::IdCollection> regionBins; ///< region idx -> snpIdxs
	std::map<uint, Utility::IdCollection> exonBinIdx; ///< region idx -> snpIdxs
	std::map<uint, Utility::IdCollection> intronBinIdx; ///< region idx -> snpIdxs
	std::map<uint, Utility::IdCollection> regulatoryBinIdx; ///< region idx -> snpIdxs
	std::map<uint, Utility::IdCollection> unknownBinIdx; ///< region idx -> snpIdxs
	std::map<uint, Utility::IdCollection> expandedRegionBins; ///< regions that are derived from groups. The collection points to the unique meta group IDs

	//The first thing we have to do is figure out
	//which SNPs are binnable and to which genes they are attributed
	/** chrom -> {block -> snpIdxs} */
	std::map<uint, std::map<uint, Utility::IdCollection> > intergenic;
	std::map<uint, std::set<uint> > regionToBinnable; ///< for any region, this is the set of indexes we are interested in
	uint snpCount = snps.Size();
	for (uint i = 0; i < snpCount; i++) {
		Utility::Locus& l = snps[i];
		if ((1.0 - l.MajorAlleleFreq()) < mafCutoff) {
			Utility::IdCollection regionIDs;
			snps.GetRegionCoverage(i, regionIDs);
			Utility::IdCollection::iterator itr = regionIDs.begin();
			Utility::IdCollection::iterator end = regionIDs.end();
			//std::cerr<<"\tBinnable Region >> "<<i<<"\t"<<snps[i].RSID()<<"\t"<<snps[i].pos<<" - \n";
			if (regionIDs.size() == 0) {
				intergenic[l.chrom][(l.pos / IntergenicBinWidth)].insert(i);
			}
			while (itr != end) {
				//std::cerr<<"\tBinnable Region : "<<i<<"\t"<<snps[i].RSID()<<"\t"<<regions[*itr].name<<"\n";
				//std::cerr<<"\t\t\t"<<*itr<<"\t"<<regions[*itr].name<<"\t"<<regions[*itr].effStart<<"-"<<regions[*itr].effEnd<<"\n";
				regionToBinnable[*itr++].insert(i);
			}
		} else {
			genotypeMap[i] = genotypeMap.size();
		}
	}
	Utility::IdCollection genesUsed;
	std::map<uint, Knowledge::GroupManagerDB>::iterator meta = groups.begin();
	std::map<uint, Knowledge::GroupManagerDB>::iterator metaEnd = groups.end();
	while (meta != metaEnd) {
		Knowledge::GroupManagerDB& gmgr = meta->second;
		std::map<uint, Utility::IdCollection> regionLookup;
		gmgr.BuildRegionCollections(regionLookup);
		CollectGroupLeaves(gmgr, regionLookup, regionToBinnable, leaves,
				genesUsed);
		meta++;
	}

	Utility::IdCollection::iterator geneNotUsed = genesUsed.end();
	//Work through all remaining genes that weren't covered by pathways and add
	//those as bins (if there are no groups, then we'll have all that have
	//some binnable SNPs)
	uint geneCount = regions.Size();
	for (uint i = 0; i < geneCount; i++) {
		if (genesUsed.find(i) == geneNotUsed)
			regionBins[i] = regionToBinnable[i];

	}

	// Break out genes into intron/exons/etc if necessary
	if (ExpandByExons) {
		std::map<uint, Utility::IdCollection>::iterator itr =
				regionBins.begin();
		std::map<uint, Utility::IdCollection>::iterator end = regionBins.end();

		while (itr != end) {
			Utility::IdCollection exons;
			Utility::IdCollection introns;
			Utility::IdCollection regulatory;
			Utility::IdCollection unknown;

			if (itr->second.size() > BinTraverseThreshold) {
				Utility::IdCollection::iterator sitr = itr->second.begin();
				Utility::IdCollection::iterator send = itr->second.end();

				while (sitr != send) {
					if (snps[*sitr].IsExon())
						exons.insert(*sitr);
					else if (snps[*sitr].IsIntron())
						introns.insert(*sitr);
					else if (snps[*sitr].IsRegulatory())
						regulatory.insert(*sitr);
					else
						unknown.insert(*sitr);

					sitr++;
				}
			}

			std::map<uint, Utility::IdCollection>::iterator prev = itr;
			itr++;

			if (exons.size() > 0 || introns.size() > 0
					|| regulatory.size() > 0) {
				if (exons.size() > 0)
					exonBinIdx[prev->first] = exons;
				if (introns.size() > 0)
					intronBinIdx[prev->first] = introns;
				if (unknown.size() > 0)
					unknownBinIdx[prev->first] = unknown;
				if (regulatory.size() > 0)
					regulatoryBinIdx[prev->first] = regulatory;
				regionBins.erase(prev);
			}
		}
	}

	// At this point, we have groups (pathways that reached the leaf point), genes 
	// and even intron/exons. These are all considered leaves and should only
	// overlap by contents at the pathway level (2 pathways might have some of the 
	// same genes present). From here, we create the bins:
	uint binID = 0;
	binIDs = std::vector<uint>(snps.Size(), (uint) -1);
	binNames.clear();
	GenerateBins(regionBins, regions, locusRemap, binID, "");
	GenerateBins(exonBinIdx, regions, locusRemap, binID, "ex");
	GenerateBins(intronBinIdx, regions, locusRemap, binID, "in");
	GenerateBins(regulatoryBinIdx, regions, locusRemap, binID, "reg");
	GenerateBins(unknownBinIdx, regions, locusRemap, binID, "unk");
	GenerateBins(leaves, regionToBinnable, locusRemap, binID);
	GenerateIntergenicBins(intergenic, locusRemap, binID);

	return std::make_pair(binID + 1, genotypeMap.size()+1);

}

void BinManager::GenerateIntergenicBins(
		std::map<uint, std::map<uint, Utility::IdCollection> >& intergenic,
		std::vector<uint>& indexLookup,
		uint& binID) {
	std::map<uint, std::map<uint, Utility::IdCollection> >::iterator chrom =
			intergenic.begin();
	std::map<uint, std::map<uint, Utility::IdCollection> >::iterator chromEnd =
			intergenic.end();

	while (chrom != chromEnd) {
		std::map<uint, Utility::IdCollection>::iterator seg =
				chrom->second.begin();
		std::map<uint, Utility::IdCollection>::iterator segEnd =
				chrom->second.end();
		std::string chromosome = Utility::ChromFromIntChr(chrom->first);
		while (seg != segEnd) {
			GenerateBin(seg->second, indexLookup, binID,
					chromosome + "-" + Utility::ToString(seg->first));
			seg++;
		}
		chrom++;
	}
}

void BinManager::GenerateBins(std::map<uint, Utility::IdCollection>& binData,
Knowledge::RegionManagerDB& regions,
std::vector<uint>& indexLookup,
uint &binID,
const std::string& type) {
	std::map<uint, Utility::IdCollection>::iterator itr = binData.begin();
	std::map<uint, Utility::IdCollection>::iterator end = binData.end();

	while (itr != end) {
		std::string binName = regions[itr->first].name;
		if (type.length() > 0)
			binName += "-" + type;
		GenerateBin(itr->second, indexLookup, binID, binName);
		itr++;
	}
}

void BinManager::GenerateBin(Utility::IdCollection& binData,
		std::vector<uint>& indexLookup, uint &binID, const std::string& name) {
	if (binData.size() >= MinBinSize) {
		Utility::IdCollection::iterator sItr = binData.begin();
		Utility::IdCollection::iterator sEnd = binData.end();

		while (sItr != sEnd) {
			uint snpIdx = *sItr++; ///< This is the index at SNP data

			if (binIDs[snpIdx] != (uint) -1) {
				multiBins.insert(std::make_pair(snpIdx, binIDs[snpIdx]));
				multiBins.insert(std::make_pair(snpIdx, binID));
			} else
				binIDs[snpIdx] = binID;

		}
		binNames.push_back(name);
		binID++;
	}

}
void BinManager::GenerateBins(
		std::map<Knowledge::Group*, Utility::IdCollection>& groups,
		std::map<uint, std::set<uint> > binnables,
		std::vector<uint>& indexLookup,
		uint &binID) {
	std::map<Knowledge::Group*, Utility::IdCollection>::iterator itr =
			groups.begin();
	std::map<Knowledge::Group*, Utility::IdCollection>::iterator end =
			groups.end();

	while (itr != end) {
		Knowledge::Group* group = itr->first;
		GenerateBin(binnables[group->id], indexLookup, binID, group->name);
		itr++;
	}
}

void BinManager::CollectGroupLeaves(Knowledge::GroupManagerDB& gmgr,
		std::map<uint, Utility::IdCollection>& regionLookup,
		std::map<uint, std::set<uint> >& regionToBinnable,
		std::map<Knowledge::Group*, Utility::IdCollection>& leaves,
		Utility::IdCollection& genesUsed,
		Utility::IdCollection& visited,
		uint groupIdx) {

	if (visited.find(groupIdx) == visited.end()) {
		visited.insert(groupIdx);
		Utility::IdCollection& geneIDs = regionLookup[groupIdx];

		Utility::IdCollection snps;
		Utility::IdCollection::iterator geneIdx = geneIDs.begin();
		Utility::IdCollection::iterator geneEnd = geneIDs.end();

		while (geneIdx != geneEnd) {
			snps.insert(regionToBinnable[*geneIdx].begin(),
					regionToBinnable[*geneIdx].end());
			geneIdx++;
		}
		Knowledge::Group& group = gmgr[groupIdx];
		// If we are still too big and we have children...
		if (snps.size() > BinTraverseThreshold && group.groups.size() > 0) {
			Utility::IdCollection::iterator groupItr = group.groups.begin();
			Utility::IdCollection::iterator groupEnd = group.groups.end();

			while (groupItr != groupEnd)
				CollectGroupLeaves(gmgr, regionLookup, regionToBinnable, leaves,
						genesUsed, visited, *groupItr++);
		} else {
			genesUsed.insert(group.regions.begin(), group.regions.end());
			leaves[&group] = snps;
		}

	}
}
void BinManager::CollectGroupLeaves(Knowledge::GroupManagerDB& gmgr,
		std::map<uint, Utility::IdCollection>& regionLookup,
		std::map<uint, std::set<uint> >& regionToBinnable,
		std::map<Knowledge::Group*, Utility::IdCollection>& leaves,
		Utility::IdCollection& genesUsed) {

	Utility::IdCollection visited;
	Utility::IdCollection root = gmgr.Root();
	Utility::IdCollection::iterator itr = root.begin();
	Utility::IdCollection::iterator end = root.end();

	while (itr != end)
		CollectGroupLeaves(gmgr, regionLookup, regionToBinnable, leaves,
				genesUsed, visited, *itr++);

}

std::set<uint> BinManager::ParseSNP(uint snpIndex, std::vector<char>& genotypes,
		std::vector<Individual>& data) {
	std::multimap<uint, uint>::iterator notMulti = multiBins.end();
	std::multimap<uint, uint>::iterator multiStart = multiBins.lower_bound(
			snpIndex);
	std::multimap<uint, uint>::iterator multiStop = multiBins.upper_bound(
			snpIndex);

	bool isGenotype = binIDs[snpIndex] == (uint) -1;
	uint gtIndex = 0;
	uint binIndex = 0;
	if (isGenotype)
		gtIndex = genotypeMap[snpIndex];
	else
		binIndex = binIDs[snpIndex];
	bool isMulti = multiBins.find(snpIndex) != notMulti;
	uint genotypeCount = genotypes.size();

	std::set<uint> bins;

	for (uint i = 0; i < genotypeCount; i++) {
		if (isGenotype) {
			data[i].SetGenotype(gtIndex, genotypes[i]);
		} else {
			if (isMulti) {
				std::multimap<uint, uint>::iterator itr = multiStart;
				while (itr != multiStop) {
					bins.insert(itr->second);
					data[i].IncrementBin(itr++->second, genotypes[i]);
				}
			} else {
				bins.insert(binIDs[snpIndex]);
				data[i].IncrementBin(binIndex, genotypes[i]);
			}
		}
	}
	return bins;
}

void BinManager::DescribeLocus(uint snpIndex, std::ostream& os,
		Knowledge::RegionManagerDB& regions, Knowledge::SnpDataset& snps) {
	Utility::Locus &snp = snps[snpIndex];
	Utility::IdCollection regionIDs;
	snps.GetRegionCoverage(snpIndex, regionIDs);
	bool isRare = (1.0 - snp.MajorAlleleFreq()) < mafCutoff;
	Utility::StringArray regionNames;
	regions.RegionNames(regionIDs, regionNames);

	Utility::StringArray bins;
	if (binIDs[snpIndex] != (uint) -1) {
		if (multiBins.find(snpIndex) == multiBins.end())
			bins.push_back(binNames[binIDs[snpIndex]]);
		else {
			std::multimap<uint, uint>::iterator itr = multiBins.lower_bound(
					snpIndex);
			std::multimap<uint, uint>::iterator end = multiBins.upper_bound(
					snpIndex);
			while (itr != end) {
				bins.push_back(binNames[itr->second]);
				itr++;
			}
		}
	}
	if (isRare)
		os << "Rare ";
	os << "Variant," << Utility::Join(regionNames, ":") << ","
			<< Utility::Join(bins, ":") << "\n";

}

/*
 void BinManager::BinSNPs(Utility::IdCollection& snpIdx,
 Knowledge::SnpDataset& snps,
 std::map<uint, uint> locusRemap,
 std::string& name,
 uint &binID) {
 if (snpIdx.size() > BinTraverseThreshold && ExpandByFunction) {
 std::map<uint, Utility::IdCollection> functionList;

 Utility::IdCollection::iterator itr			= snpIdx.begin();
 Utility::IdCollection::iterator end			= snpIdx.end();
 while (itr != end) {
 functionList[snps[*itr]->Function()]	= *itr;
 itr++;
 }

 if (functionList.size() > 1) {
 std::map<uint, Utility::IdCollection>::iterator fitr = functionList.begin();
 std::map<uint, Utility::IdCollection>::iterator fend = functionList.end();
 while (fitr!= fend) {
 Utility::IdCollection::iterator fsitr = fitr->begin();
 Utility::IdCollection::iterator fsend = fitr->end();
 if (fsitr != fsend) {
 while (fsitr != fsend) {
 uint realIndex = locusRemap[*fsitr++];
 if (binIDs[realIndex] != (uint)-1) {
 multiBins[realIndex].insert(binIDs[realIndex]);
 multiBins[realIndex].insert(binID);
 } else {
 binIDs[realIndex] = binID;
 }
 }
 binNames[binID] = name + "-" + Utility::Locus::FunctionID(fitr->first);
 binID++;
 }
 fitr++;
 }
 }
 }

 */
}
