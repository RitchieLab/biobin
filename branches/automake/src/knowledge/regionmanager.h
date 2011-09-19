/* 
 * File:   regionmanager.h
 * Author: torstees
 *
 * Basic interface for loading and dereferencing region entities. This must be
 * untied to SOCI, since it might be included with programs that will only be
 * getting data from the output files, and not from the database itself
 *
 * Created on March 7, 2011, 2:33 PM
 */

#ifndef REGIONMANAGER_H
#define	REGIONMANAGER_H

#include <vector>
#include <fstream>
#include <map>
#include <utility>					// make_pair
#include "region.h"
#include "utility/exception.h"
#include "def.h"
#include "snpdataset.h"
//#include "genegenemodel.h"

namespace Knowledge {

class RegionManager {
public:
	typedef std::vector<Region> RegionCollection;
	RegionManager() {} 
	RegionManager(const RegionManager& orig) : regions(orig.regions) { }
	virtual ~RegionManager() {}

	void WriteArchive(std::ostream& file, const char *sep = "\t");

	/**
	 * Writes regions to archive (will switch to binary based on BinaryArchive)
    * @param filename
    */
	void WriteArchive(const char *filename, const char *sep = "\t");

	/**
	 * Loads regions from archive (will switch to binary, based on BinaryArchive)
    * @param filename
    */
	void LoadArchive(const char *filename, const char *sep = "\t");

	/**
	 * Adds region and returns
    * @param name		the primary name
    * @param id		the gene_id from the database
    * @param start	assigned to both true and effective start
    * @param stop		assigned to both true and effective stop
	 * @param aliases Comma separated list of aliases
    * @return returns a reference to the actual region
    */
	Region& AddRegion(const char *, uint id, char chrom, uint start, uint stop, const char *aliases = "");

	void AddMetaID(uint id, MetaGroup::Type groupType, Utility::IdCollection& regionIDs);
	void AddMetaID(uint id, MetaGroup::Type groupType, uint regionId);
	
	/**
	 * Adds region and returns it's index
    * @param name		primary name
    * @param id		gene_id
    * @param effStart	effective Start
    * @param effStop		effective stop
    * @param trueStart	true start
    * @param trueStop	true stop
	 * @param aliases Comma separated list of aliases
    * @return returns a reference to the actual region
    */
	Region& AddRegion(const char *name, uint id, char chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const char* aliases = "");

	
	void RegionNames(Utility::IdCollection& regionIDs, Utility::StringArray& regionNames);
	/**
	 * Grant access to member items using the index
    * @param idx
    * @return
    */
	Region& operator[](uint idx);

	const Region& operator[](uint idx) const;

	Region& operator[](const std::string& alias);
	//const Region& operator[](const char *alias) const;

	//Use the (id) operator to query by ID rather than index
	uint operator()(const std::string& alias);
	uint operator()(uint id);

	/**
	 * Build out the segments for a given chromosome. These are effectively the bins associated with
	 * the regions. We are assuming that each segment will represent the highest resolution. So
	 * if that needs to be variable, you'll have to add some variable to control it's 
	 * behavior.
    * @param chromosome  -- Indicates which chromosome we are looking at
    * @param segments
    */
	void BuildRegionSegments(char chromosome, RegionContainer& segments);

	void GenerateGeneReport(const char *filename, SnpDataset& dataset);
	/**
	 * Generates a lookup map snpIdx -> regionIdx
    * @param geneLookup
    */
	void BuildSnpGeneMap(std::multimap<uint, uint>& geneLookup);

	void AssociateSNPs(SnpDataset& snps);


	/**
	 * Removes any region that has no SNPs associated with it
    * @return Number of regions removed
    */
	uint Squeeze();

	/**
	 * Returns the number of regions contained within
    */
	uint Size();

	bool ValidGeneGene(uint r, uint l);

	/**
	 * Builds associations for local genes as found in dataset, snps. 
    * @param associations	snp idx->gene idx
    * @param snps
    */
	void BuildSnpAssociations(std::multimap<uint, uint>& associations, const SnpDataset& snps);
		
	/**
	 * Checks the ModelGenerationMode to determine if there is a possilibity of generating models
    * @param groupContents List of region IDs to consult
    * @return
    */
	bool DoGenerateModels(Utility::IdCollection& groupContents);
	
	static ModelGenerationMode::Type modelGenerationType;
protected:
	void WriteArchiveBinary(const char *filename);
	void LoadArchiveBinary(const char *filename);
	std::map<uint, uint> idToIndex;						///< id=>index lookup
	std::map<std::string, uint> aliasToIndex;			///< alias=>index lookup
	RegionCollection regions;
};


inline
uint RegionManager::Size() {
	return regions.size();
}

inline
bool RegionManager::ValidGeneGene(uint r, uint l) {
	return modelGenerationType != ModelGenerationMode::DD_ONLY || (regions[r].CountDDCapable() + regions[l].CountDDCapable() > (uint)0);
}
inline
bool RegionManager::DoGenerateModels(Utility::IdCollection& groupContents) {
	if (modelGenerationType == ModelGenerationMode::ALL_MODELS)
		return true;
	Utility::IdCollection::iterator itr = groupContents.begin();
	Utility::IdCollection::iterator end = groupContents.end();

	uint ddCapable = 0;
	//Build up the group collection data
	while (itr != end) 
		ddCapable += regions[*itr++].CountDDCapable();

	return ddCapable > 0;
}

/**
 * Functions similarly to the regular functionality, except it returns
 * the associations instead of building them into the system. We do this
 * for things like getting a list of genes associated with a speclialized
 * set of SNPs. In this situation, we'll only get hits from genes that
 * are observed in the real dataset...
 * @param associations
 */
inline
void RegionManager::BuildSnpAssociations(std::multimap<uint, uint>& associations, const SnpDataset& snps) {
	Utility::IdCollection snplist;
	int count = regions.size();
	for (int i=0; i<count; i++) {
		Region& region = regions[i];
		snplist.clear();
		snps.RangeSnpLookup(region.chrom, region.effStart, region.effEnd, snplist);
		Utility::IdCollection::iterator sitr = snplist.begin();
		Utility::IdCollection::iterator send = snplist.end();
		
		while (sitr != send) {
			associations.insert(std::make_pair(*sitr, i));
			sitr++;
		}
	}
}
inline
void RegionManager::BuildRegionSegments(char chromosome, RegionContainer& segments) {
	RegionCollection::iterator itr = regions.begin();
	RegionCollection::iterator end = regions.end();

	uint count = 0;
	uint i=0;								///< Gene indexes
	while (itr != end) {
		Knowledge::Region &r =*itr++;
		//std::cerr<<r.name<<"\t"<<(uint)r.chrom<<"\n";
		if (r.chrom == chromosome) {
			count++;
			segments.AddSegment(r.effStart, r.effEnd, i);
		}
		i++;
	}
	
}


inline
void RegionManager::RegionNames(Utility::IdCollection& ids, Utility::StringArray& regionNames) {
	Utility::IdCollection::iterator ritr = ids.begin();
	Utility::IdCollection::iterator rend = ids.end();
		while (ritr != rend) 
			regionNames.push_back(regions[*ritr++].name);
}

inline
void RegionManager::GenerateGeneReport(const char *filename, SnpDataset& dataset) {
	std::ofstream file(filename);
	RegionCollection::iterator itr = regions.begin();
	RegionCollection::iterator end = regions.end();
	file<<"Gene Name,Eff. Start,Eff. Stop,Alias List,Start,End";
	if (Knowledge::SnpDataset::detailedReport)
		file<<",SNPs\n";
	while (itr != end) {
		Knowledge::Region &r =*itr++;
		file<<r.name<<","
			 <<Utility::ChromFromInt(r.chrom - 1)<<","
			 <<r.effStart<<","
			 <<r.effEnd<<","
			 <<r.trueStart<<","
			 <<r.trueEnd<<","
			 <<r.GetAliasString(":");
		if (Knowledge::SnpDataset::detailedReport) {
			file<<","<<r.GetSnpString(":", dataset);
		}
		file<<"\n";
	}
}

inline
uint RegionManager::operator()(uint id) {
	if (idToIndex.find(id) != idToIndex.end())
		return idToIndex[id];
	return (uint)-1;
}

inline
uint RegionManager::operator()(const std::string& alias)  {
	if (aliasToIndex.find(alias) != aliasToIndex.end())
		return aliasToIndex[alias];
	return (uint)-1;
}

inline
void RegionManager::AddMetaID(uint id, MetaGroup::Type groupType, Utility::IdCollection& ids) {
	Utility::IdCollection::iterator itr = ids.begin();
	Utility::IdCollection::iterator end = ids.end();

	while (itr != end) {
		regions[*itr++].AddMetaID(groupType, id);
	}
}

inline
void RegionManager::AddMetaID(uint id, MetaGroup::Type groupType, uint regionId) {
	regions[regionId].AddMetaID(groupType, id);
}

inline
void RegionManager::WriteArchive(std::ostream& file, const char *sep) {
	file<<"Region Name"<<sep<<"True Begin"<<sep<<"True End"<<sep<<"Eff. Begin"<<sep<<"Eff. End"<<sep<<"Groups"<<sep<<"Aliases"<<sep<<"SNPs"<<"\n";

	for (RegionCollection::iterator itr=regions.begin(); itr<regions.end(); itr++) {
		itr->WriteToArchive(file, sep);
	}

}
inline
void RegionManager::WriteArchive(const char* filename, const char *sep) {
	if (BinaryArchive)
		WriteArchiveBinary(filename);
	else {
		std::ofstream file(filename);
		WriteArchive(file, sep);
	}
}

inline
uint RegionManager::Squeeze() {
	uint count = regions.size();
	Utility::IdCollection empties;

	for (uint i=0; i<count; i++) {
		if (regions[i].SnpCount() == 0)
			empties.insert(i);
	}

	RegionCollection backup = regions;
	regions.clear();
	std::map<uint, uint> remap;

	for (uint i=0; i<count; i++) {
		if (empties.find(i) == empties.end()) {
			remap[i] = regions.size();
			regions.push_back(backup[i]);
		} else {
			remap[i] = (uint)-1;
		}
	}

	std::map<std::string, uint>::iterator itr = aliasToIndex.begin();
	std::map<std::string, uint>::iterator end = aliasToIndex.end();
	while (itr != end) {
		itr->second = remap[itr->second];
		itr++;
	}

	std::map<uint, uint>::iterator iitr = idToIndex.begin();
	std::map<uint, uint>::iterator iend = idToIndex.end();
	while (iitr != iend) {
		iitr->second = remap[iitr->second];
		iitr++;
	}
	return empties.size();
}

inline
void RegionManager::WriteArchiveBinary(const char *filename) {
	std::ofstream file(filename, std::ios::binary);
	uint count = regions.size();
	file.write((char*)&count, 4);


	for (RegionCollection::iterator itr=regions.begin(); itr<regions.end(); itr++) {
		itr->WriteToArchiveBinary(file);
	}
}

inline
void RegionManager::LoadArchive(const char *filename, const char *sep) {
	if (BinaryArchive)
		LoadArchiveBinary(filename);
	else {
		std::ifstream file(filename);
		char line[4096];
		file.getline(line, 4096);

		Region reg;
		while (file.good() && !file.eof()) {
			if (reg.LoadFromArchive(file, sep))
				regions.push_back(reg);
		}
	}
}

inline
void RegionManager::LoadArchiveBinary(const char *filename) {
	std::ifstream file(filename, std::ios::binary);

	uint count = 0;
	file.read((char*)&count, 4);

	Region reg;
	for (uint i=0; i<count; i++) {
		reg.LoadFromArchiveBinary(file);
		regions.push_back(reg);
	}
}

inline
void RegionManager::AssociateSNPs(SnpDataset& snps) {
	RegionCollection::iterator itr = regions.begin();
	RegionCollection::iterator end = regions.end();

	Utility::IdCollection snplist;
	while (itr != end) {
		snplist.clear();
		std::string geneName = itr->name;

		assert(itr->chrom > 0);
		snps.RangeSnpLookup(itr->chrom, itr->effStart, itr->effEnd, snplist);
		itr++->AddSNPs(snplist);
	}
	
	if (snps.Size() > 0)	
		Squeeze();
	
	itr = regions.begin();
	end = regions.end();
	uint i = 0;
	while (itr != end) {
		snps.AddRegion(itr->chrom, itr->effStart, itr->effEnd, i++);
		itr++;
	}
}

inline
Region& RegionManager::AddRegion(const char *name, uint id, char chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const char* aliases) {
	uint idx = regions.size();
	assert(chrom > 0);
	Region reg(name, id, chrom, effStart, effStop, trueStart, trueStop);
	reg.AddAliases(aliases);
	regions.push_back(reg);
	idToIndex[id] = idx;
	aliasToIndex[name] = idx;
	
	Utility::StringArray aliasList = Utility::Split(aliases, ",");
	Utility::StringArray::iterator itr = aliasList.begin();
	Utility::StringArray::iterator end = aliasList.end();
	while (itr != end) 
		aliasToIndex[*itr++] = idx;

	return regions[idx];
}

inline
Region& RegionManager::AddRegion(const char *name, uint id, char chrom, uint start, uint stop, const char *aliases) {
	return AddRegion(name, id, chrom, start, stop, start, stop, aliases);
}

inline
Region& RegionManager::operator[](const std::string& alias) {
	if (aliasToIndex.find(alias) != aliasToIndex.end())
		return regions[aliasToIndex[alias]];
	throw Utility::Exception::IndexOutOfRange("regions", aliasToIndex[alias], regions.size());
}

inline
void RegionManager::BuildSnpGeneMap(std::multimap<uint, uint>& geneLookup) {
	uint count = regions.size();

	for (uint i=0; i<count; i++) {
		Utility::IdCollection& snpIDs = regions[i].snps;
		Utility::IdCollection::iterator itr = snpIDs.begin();
		Utility::IdCollection::iterator end = snpIDs.end();

		while (itr != end)
			geneLookup.insert(std::pair<uint,uint>(*itr++, i));
	}
}


inline
Region& RegionManager::operator[](uint idx) {
	return regions[idx];
}

inline
const Region& RegionManager::operator[](uint idx) const {
	return regions[idx];
}

}

#endif	/* REGIONMANAGER_H */

