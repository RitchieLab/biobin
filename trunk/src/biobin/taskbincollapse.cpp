/* 
 * File:   taskbincollapse.cpp
 * Author: torstees
 * 
 * Created on July 5, 2011, 12:12 PM
 */

#include "taskbincollapse.h"
#include "taskfilegeneration.h"		// for the common delimiter
namespace BioBin {
namespace Task {


/**
 This is an arbitrary choice by Eric Torstenson
 */
int	BinCollapse::maxSnpCount					= 200;
bool	BinCollapse::VisualizeGroupTrees			= true;
bool	BinCollapse::WriteKnowledgeBins			= true;

void BinCollapse::CollectSNPs(std::map<uint, std::set<uint> >& snpCollection, uint idx, Knowledge::GroupManager& mgr, Knowledge::RegionManager& regions) {
	if (snpCollection.find(idx) == snpCollection.end()) {
		Knowledge::Group& g									= mgr[idx];
		std::multimap<uint, uint>::iterator missing	= binContents.end();
		Utility::IdCollection::iterator itr				= g.regions.begin();
		Utility::IdCollection::iterator end				= g.regions.end();
		std::set<uint> localSNPs;
		// Build up the snp counts for the local object
		while (itr != end) {
			if (binContents.find(*itr)!= missing) {
				snps+=binContents.count(*itr);
				std::multimap<uint, uint>::iterator first = binContents.lower_bound(*itr);
				std::multimap<uint, uint>::iterator last = binContents.upper_bound(*itr);
				while (first != last)  {
					localSNPs.insert(first++->second);
				}
			}
			itr++;
		}
		
		snpCollection[idx] = localSNPs;
		
		// Add in counts for any children as well
		Utility::IdCollection::iterator gitr			= g.groups.begin();
		Utility::IdCollection::iterator gend			= g.groups.end();
		std::map<uint, std::set<uint> >::iterator groupMissing  = snpCollection.end();
		while (gitr != gend) {
			if (*gitr != idx) {
				if (snpCollection.find(*gitr) == groupMissing)
					CollectSNPs(snpCollection, *gitr, mgr, regions);
				std::set<uint> &snps = snpCollection[*gitr];
				localSNPs.insert(snps.begin(), snps.end());
			}
			gitr++;
		}
		snpCollection[idx] = localSNPs;
	}

}


void BinCollapse::EvaluateGroup(uint tabCount, std::map<uint, std::set<uint> >& snpCollection, std::set<uint>& visited, uint idx, Knowledge::GroupManager& mgr, Knowledge::RegionManager& regions) {
	// Iterate over each region in the group and add all of the SNPs from the binContents object
	// if it's present.   
	std::stringstream ss;
	if (visited.find(idx) != visited.end())
		return;
	
	visited.insert(idx);
	uint snpCount = snpCollection.find(idx)->second.size();
	if (snpCount == 0)
		return;
	
	Knowledge::Group &group = mgr[idx];
	if (VisualizeGroupTrees)
		visualization<<std::string(tabCount, '\t')<<group.name<<":"<<group.id<<" ("<<snpCount<<"):\n";
	// Then, if number of SNPs is small enough, write a leaf style line. (or if we don't have any children...
	
	// Right now, we have a problem with circular references. If A includes B which includes A
	// neither is successfully resolved as a leaf because group.groups.size() is greater than 0
	bool isLeaf = maxSnpCount == -1 || snpCount <= (uint)maxSnpCount;

	std::multimap<uint, uint>::iterator missing = binContents.end();

	Utility::IdCollection::iterator itr = group.regions.begin();
	Utility::IdCollection::iterator end = group.regions.end();
	if (isLeaf || group.groups.size() == 0) {
		Utility::IdCollection binIDs;
		
		Utility::IdCollection otherIDs;				///< This is for debugging 
		//std::cerr<<"ASDF: binIndex: "<<Utility::MapJoin(binIndex)<<"\n";
		while (itr != end) {
		//for (uint i=0; i<count; i++) {
			if (binContents.find(*itr)!= missing) {
				std::multimap<uint, uint>::iterator first = binContents.lower_bound(*itr);
				std::multimap<uint, uint>::iterator last = binContents.upper_bound(*itr);
				//std::cerr<<">>  ASDF: "<<group.name<<" : "<<first->first<<" => "<<first->second<<" ("<<binIndex[first->first]<<")\n";
				binIDs.insert(binIndex[first->first]);
				otherIDs.insert(first->first);
				visualization<<std::string(tabCount+1, '\t')<<" "<<regions[*itr].name<<"("<<binContents.count(*itr)<<"):\n";
				while (first != last)  {
					if (VisualizeGroupTrees) 
						visualization<<std::string(tabCount+2, '\t')<<app->Locus(first++->second).RSID()<<"\n";
				}
			}
			itr++;
		}
		if (isLeaf) {
			groupParticipant.push_back(std::make_pair(&group, binIDs));
			//std::cerr<<"ASDFASDF Group => IDs:: "<<group.name<<" :: "<<Utility::Join(binIDs, "\t")<<"\n";
			mergedBins.insert(binIDs.begin(), binIDs.end());
		}
	}
	else {
	// otherwise, write a summary and call this on the children
		uint count = group.regions.size();
		count = group.groups.size();
		Utility::IdCollection::iterator itr = group.groups.begin();
		Utility::IdCollection::iterator end = group.groups.end();
		while (itr != end) 
			EvaluateGroup(tabCount+1, snpCollection, visited, *itr++, mgr, regions);
	}
}


void BinCollapse::WriteBinData() {
	const char *sep = GenerateFiles::OutputDelimeter.c_str();
	std::string filename = app->AddReport("data", "kbins", "Knowledge Based Bin Data");
	ofstream binfile(filename.c_str());
	
	uint newBinCount = 0;
	std::map<uint, uint> binConverter;			///< Used to reposition the bin data to it's new home

	
	const std::vector<Individual> &individuals = app->Individuals();
	std::vector<Individual>::const_iterator itr = individuals.begin();
	std::vector<Individual>::const_iterator end = individuals.end();

	uint statusCount = 0;
	
	{
		std::set<float> uniqueStatus;
		while (itr != end) 
			uniqueStatus.insert(itr++->status);
		statusCount = uniqueStatus.size();
	}
	itr														= individuals.begin();
	
	
	{	// I'm going to let the names and maxhits fall out of scope...
		std::vector<uint> maxHits;				///< This is the actual max bin hits including the knowledge based bins
		std::vector<uint> origBinHits;		///< These are the original bins
		app->GetMaxBinHits(origBinHits);
		Utility::StringArray names;
		names.reserve(groupParticipant.size());
		std::vector<std::pair<Knowledge::Group*, Utility::IdCollection> >::iterator gitr = groupParticipant.begin();
		std::vector<std::pair<Knowledge::Group*, Utility::IdCollection> >::iterator gend = groupParticipant.end();
		//std::cerr<<"ASDFASDF: "<<groupParticipant.size()<<"\n";
		while (gitr != gend) {
			Knowledge::Group *g = gitr->first;
			Utility::IdCollection::iterator itr = gitr->second.begin();
			Utility::IdCollection::iterator end = gitr->second.end();
			
			uint data = 0;
			while (itr != end )
				data+=(origBinHits[*itr++]*2);
			maxHits.push_back(data);
			names.push_back(g->name);
			gitr++;
		}

		{
			Utility::StringArray binNames = app->BinNames();
			uint count = binNames.size()+1;
			Utility::IdCollection::iterator notMerged = mergedBins.end();			
			
			//std::cerr<<"ASDFASDF: Merged Bin Indexes: "<<Utility::Join(mergedBins, " ")<<"\n";
			
			for (uint i=1; i<count;i++) {
				if (mergedBins.find(i) == notMerged) {
					maxHits.push_back(origBinHits[binIndex[i]]);
					binConverter[i]			= names.size();
					names.push_back(binNames[i-1]);
				} 
			}
		}
		binfile<<"ID"<<sep<<"Status"<<sep<<"Intergenic"<<sep<<Utility::Join(names, sep)<<"\n";
		binfile<<"Totals"<<sep<<statusCount<<sep<<origBinHits[0]*2<<sep<<Utility::Join(maxHits, sep)<<"\n";
		
		

		newBinCount = maxHits.size();
	}
	
	

	std::vector<uint> binData;
	binData.reserve(newBinCount);
	while (itr != end) {
		binData.clear();
		std::vector<std::pair<Knowledge::Group*, Utility::IdCollection> >::iterator gitr = groupParticipant.begin();
		std::vector<std::pair<Knowledge::Group*, Utility::IdCollection> >::iterator gend = groupParticipant.end();
		while (gitr != gend) {
			uint data = 0;
			Utility::IdCollection::iterator bitr = gitr->second.begin();
			Utility::IdCollection::iterator bend = gitr->second.end();
			while (bitr != bend ) 
				data+=itr->BinCount(*(bitr++));
			binData.push_back(data);
			gitr++;	
		}
		std::map<uint, uint>::iterator bitr = binConverter.begin();
		std::map<uint, uint>::iterator bend = binConverter.end();
		while (bitr != bend) {
			assert(bitr->second==binData.size());
			binData.push_back(itr->BinCount(bitr->first));
			bitr++;
		}

		std::vector<uint> binHits;
		app->GetMaxBinHits(binHits);					///< This is number of SNPs
		std::vector<uint>::iterator bhitr			= binHits.begin();
		std::vector<uint>::iterator bhend			= binHits.end();
		uint count = 0;
		while (bhitr!=bhend) {
			count += *bhitr;
			bhitr++;
		}
		
		
		binfile<<itr->indID<<sep<<itr->status<<sep<<count<<sep<<Utility::Join(binData, sep)<<"\n";
		itr++;
	}

}

void BinCollapse::ExecuteTask() {
	std::cerr<<"Executing Bin Collapse\n";
	std::vector<uint> mgrKeys = app->ManagerIDs();
	std::vector<uint>::iterator itr = mgrKeys.begin();
	std::vector<uint>::iterator end = mgrKeys.end();
	binIndex							= app->GetBinLookup();

	
	((BinApplication*)app)->GenerateBinContentLookup(binContents);
	while (itr != end) {
		Knowledge::GroupManagerDB& mgr = app->GroupManager(*itr++);
		std::stringstream ss;
		ss<<"collapsed-bin-report-"<<mgr.name;
		std::string filename  = app->AddReport(ss.str().c_str(), "txt", "Hierarchical bin representation");
		if (VisualizeGroupTrees) {
			visualization.close();
			visualization.open(filename.c_str());
		}
		uint count = mgr.Size();
		std::map<uint, std::set<uint> > snpCollection;
		std::set<uint> visited;
		std::set<uint> totalSNPs;
		for (uint g=0; g<count; g++) {
			CollectSNPs(snpCollection, g, mgr, *(app->GetRegions()));
			
			std::set<uint> snps = snpCollection[g];
			totalSNPs.insert(snps.begin(), snps.end());

		}
		if (VisualizeGroupTrees)
			visualization<<mgr.name<<"\t("<<totalSNPs.size()<<"):\n";
		for (uint g=0; g<count; g++) {
			EvaluateGroup(1, snpCollection, visited, g, mgr, *(app->GetRegions()));
		}
	}
	
	visualization.close();
	//Now let's generate the binfile itself
	if (WriteKnowledgeBins) {
		WriteBinData();
	}	
	//This is going to be pretty big...so, let's not keep it around
	binContents.clear();
	binIndex.clear();
	groupParticipant.clear();
	mergedBins.clear();
}


}
}

