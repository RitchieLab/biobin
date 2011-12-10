/* 
 * File:   binmanager.cpp
 * Author: torstees
 * 
 * Created on July 21, 2011, 3:54 PM
 */

#include "binmanager.h"

#include "Bin.h"

#include "knowledge/RegionCollection.h"
#include "knowledge/Region.h"

using Knowledge::RegionCollection;
using Knowledge::Region;
using std::make_pair;

namespace BioBin {

uint BinManager::IntergenicBinWidth = 50000;
uint BinManager::BinTraverseThreshold = 50;
uint BinManager::MinBinSize = 1;
bool BinManager::ExpandByExons = true;
bool BinManager::ExpandByFunction = true;
float BinManager::mafCutoff = 0.05;
uint BinManager::maxSnpCount = 200;

BinManager::BinManager() : _total_variants(0){
}

BinManager::~BinManager() {
	set<Bin*>::const_iterator itr = _bin_list.begin();
	set<Bin*>::const_iterator end = _bin_list.end();
	while(itr != end){
		delete *itr;
		++itr;
	}
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
void BinManager::InitBins(
		const map<uint, Knowledge::GroupCollection*> &groups,
		const Knowledge::RegionCollection& regions,
		const vector<Knowledge::Locus*>& loci) {
	/**
	 * The following items will result in a bin once we reach the end of the 
	 * function. So, if we refine an item at one level, we should remove it from
	 * the original structure. These will be the index into the dataset-not into
	 * the vcf file (the locusRemap is used at the very last step only)
	 */
	/** Group -> snpIdx*/

	vector<Knowledge::Locus*>::const_iterator l_itr = loci.begin();
	vector<Knowledge::Locus*>::const_iterator l_end = loci.end();

	_total_variants = loci.size();

	while(l_itr != l_end){
		Knowledge::Locus& l = **l_itr;
		if ((1.0 - l.majorAlleleFreq()) < mafCutoff) {

			_rare_variants.insert(&l);

			// First, find all of the regions that contain this locus
			RegionCollection::const_region_iterator r_itr =
					regions.positionBegin(l.getChrom(), l.getPos());
			RegionCollection::const_region_iterator r_end =
					regions.positionEnd(l.getChrom(), l.getPos());

			Bin* curr_bin;
			if (r_itr == r_end){
				//Add to intergenic
				pair<short, int> key = make_pair(l.getChrom(), l.getPos() / IntergenicBinWidth);
				map<pair<short, int>, Bin*>::const_iterator i_bin =
						_intergenic_bins.find(key);

				if (i_bin == _intergenic_bins.end()){
					// OK, this bin is nonexistent
					curr_bin = new Bin(key.first, key.second);
					_intergenic_bins.insert(make_pair(key, curr_bin));
					_bin_list.insert(curr_bin);
				}else{
					curr_bin = (*i_bin).second;
				}
				curr_bin->addLocus(&l);
				_locus_bins[&l].insert(curr_bin);
			}
			while (r_itr != r_end){
				// For each region, get all of the groups it is in
				Region::const_group_iterator g_itr = (*r_itr)->groupBegin();
				Region::const_group_iterator g_end = (*r_itr)->groupEnd();

				// If not in any groups, add to region bins
				if(g_itr == g_end){
					int id = (*r_itr)->getID();
					map<int, Bin*>::const_iterator rm_itr = _region_bins.find(id);
					if (rm_itr == _region_bins.end()){
						curr_bin = new Bin(*r_itr);
						_bin_list.insert(curr_bin);
						_region_bins.insert(make_pair(id, curr_bin));
					}else{
						curr_bin = (*rm_itr).second;
					}
					curr_bin->addLocus(&l);
					_locus_bins[&l].insert(curr_bin);
				}

				// add to all group bins that it is a member of
				while(g_itr != g_end){
					int id = (*g_itr)->getID();
					map<int, Bin*>::const_iterator gm_itr = _group_bins.find(id);
					if (gm_itr == _group_bins.end()){
						curr_bin = new Bin(*g_itr);
						_bin_list.insert(curr_bin);
						_group_bins.insert(make_pair(id,curr_bin));
					}else{
						curr_bin = (*gm_itr).second;
					}
					curr_bin->addLocus(&l);
					_locus_bins[&l].insert(curr_bin);
					++g_itr;
				}
				++r_itr;
			}
		}
		++l_itr;
	}

	// At this point, we have all of the top level bins constructed and
	// stored in the variable _bin_list.  We should now collapse the
	// bins according to the preferences given
	collapseBins();
}

void BinManager::printBins(std::ostream& os, Knowledge::Locus* l,
		const string& sep){
	map<Knowledge::Locus*, set<Bin*> >::const_iterator m_itr = _locus_bins.find(l);

	if(m_itr != _locus_bins.end()){
		set<Bin*>::const_iterator s_itr = (*m_itr).second.begin();
		set<Bin*>::const_iterator s_end = (*m_itr).second.end();
		if (s_itr != s_end){
			os << (*s_itr)->getName();
			while(++s_itr != s_end){
				os << sep << (*s_itr)->getName();
			}
		}
	}

}

void BinManager::collapseBins(){
	return;
}

} // namespace BioBin;




/*	for (uint i = 0; i < snpCount; i++) {

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

 */


/*void BinManager::GenerateIntergenicBins(
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
		std::string chromosome = Utility::ChromFromIntChr(chrom->first-1);
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
 */
/* THIS WAS COMMENTED OUT PREVIOUSLY!
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
