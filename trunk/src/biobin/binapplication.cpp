/* 
 * File:   binapplication.cpp
 * Author: torstees
 * 
 * Created on June 22, 2011, 10:35 AM
 */

#include "binapplication.h"

#include <iostream>

#include "knowledge/Locus.h"
#include "knowledge/Region.h"
#include "knowledge/Group.h"

#include "knowledge/liftover/ConverterSQLite.h"

namespace BioBin {

BinApplication::BinApplication() {
}


BinApplication::~BinApplication() {
}



void BinApplication::InitBins() {
	
	//Utility::IdCollection variants;
	//Utility::IdCollection rareVariants;

	//if (vcfimporter.Open(filename.c_str(), (char)-1)) {
		//We now have our snp dataset set up-so it's time to start the binning process
	//std::pair<uint, uint> binGenotypeCounts;
	//binGenotypeCounts = binData.InitBins(groups, regions, dataset);

	binData.InitBins(groups, *regions, dataset);

	//binData.CollectVariantGroups(variants, rareVariants);
	std::cerr<<"   Total SNPS:   "<<std::setw(10)<<std::right<<dataset.size()<<"\n"
				<<"   Variants:     "<<std::setw(10)<<std::right<<binData.numVariants()<<"\n"
				<<" * Rare Variants:"<<std::setw(10)<<std::right<<binData.numRareVariants()<<"\n"
				<<"   Total Bins:   "<<std::setw(10)<<std::right<<binData.numBins()<<"\n";

	std::cerr<<"\n   * Rare variants are those whose minor alleles sum is below: "<<BinManager::mafCutoff<<"\n";

	if (binData.numBins() < 500) {
		std::cerr<<"\n\nBin Name\tSNP Count\n";
		BinManager::const_iterator itr = binData.begin();
		BinManager::const_iterator end = binData.end();
		while(itr != end){
			std::cerr << (*itr)->getName() << "\t" << (*itr)->getSize() << "\n";
		}
	}
	/*
	Utility::IdCollection::iterator varItr = variants.begin();
	Utility::IdCollection::iterator varEnd = variants.end();
	while (varItr != varEnd) 
		GenotypeStorage::alleleCount.push_back(dataset[*varItr++].alleles.size());


	Utility::StringArray individualIDs = vcfimporter.GetIndividualIDs();
	Utility::StringArray::iterator iitr = individualIDs.begin();
	Utility::StringArray::iterator iend = individualIDs.end();

	uint individualCount = individualIDs.size();
	individuals = std::vector<Individual>(individualCount);

	uint i=0;
	while (iitr != iend) {
		individuals[i++].Init(*iitr, binGenotypeCounts.second, binGenotypeCounts.first + 1);
		iitr++;
	}		
*/
	std::string ofn = AddReport("locus", "csv", "Locus Description");
	std::ofstream locusFile(ofn.c_str());
	locusFile<<"Chromosome,Location,ID,type,gene(s),bin name(s)\n";
	vector<Knowledge::Locus*>::const_iterator itr = dataset.begin();
	vector<Knowledge::Locus*>::const_iterator end = dataset.end();
	while(itr != end){
		locusFile << (**itr) << ",";
		locusFile << string(((1.0 - (*itr)->majorAlleleFreq()) < BinManager::mafCutoff) ? "Rare " : "");
		locusFile << string("Variant") << string(",");
		// Print the genes here
		Knowledge::RegionCollection::const_region_iterator r_itr =
				regions->positionBegin((*itr)->getChrom(), (*itr)->getPos());
		Knowledge::RegionCollection::const_region_iterator r_end =
				regions->positionBegin((*itr)->getChrom(), (*itr)->getPos());
		if (r_itr != r_end){
			locusFile << (*itr)->getID();
			while(++r_itr != r_end){
				locusFile << ":" << (*itr)->getID();
			}
		}
		locusFile << ",";
		// Now print the bins
		binData.printBins(locusFile, *itr);
		locusFile << "\n";
		++itr;
	}
	locusFile.close();

	/*
	//binIDs.resize(locusArray.size(), (uint)-1);
	for (uint i=0; i<locusCount; i++) {
		std::vector<char> genotypes(individualCount, (char)-1);
		if (dataset[i].chrom > 0) {
			vcfimporter.ParseSNP(locusRemap[i], genotypes);
			binData.ParseSNP(i, genotypes, individuals);
		}
	}

	//We should have binned data and genotypes sorted out
	ApplyPhenotypes();

	//}
	return binGenotypeCounts;
	*/
}

void BinApplication::writeBinData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	_pop_mgr.printBins(file, binData, sep);
	file.close();
}

void BinApplication::writeGenotypeData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	_pop_mgr.printGenotypes(file, sep);
	file.close();
}


/*
void BinApplication::GetMaxBinHits(std::vector<uint>& hits) {
	std::vector<std::set<uint> > binContributors;
	binData.BuildContributorList(binContributors);
	uint binCount = binContributors.size();
	hits = std::vector<uint>(binCount, 0);

	for (uint i=0; i<binCount; i++) {
		hits[i] = binContributors[i].size();
	}
}

Utility::Locus& BinApplication::Locus(uint idx) {
	return dataset[idx];
}

std::map<uint, uint> BinApplication::GetBinLookup() {
	return binIndex;
}

const std::vector<Individual>& BinApplication::Individuals() {
	return individuals;
}

const Utility::StringArray& BinApplication::BinNames() {
	return binData.BinNames();
}

const Knowledge::Region& BinApplication::GetRegion(uint idx) {
	return regions[idx];
}

void BinApplication::GetBinContributors(std::vector<std::set<uint> >& contributors) {
	binData.BuildContributorList(contributors);
}


*/
/*
void BinApplication::ApplyPhenotypes() {
	Utility::StringArray::iterator itr = phenotypeFilenames.begin();
	Utility::StringArray::iterator end = phenotypeFilenames.end();
	std::map<std::string, std::string> phenotypeLookup;
	while (itr != end) {
		Utility::StringArray ids;
		std::string contents = Utility::LoadContents(itr->c_str());
		ids = Utility::Split(contents.c_str(), "\n");
		Utility::StringArray::iterator id = ids.begin();
		Utility::StringArray::iterator kvend = ids.end();

		while (id != kvend) {
			Utility::StringArray kv = Utility::Split(id->c_str());
			if (kv.size() > 1)
				phenotypeLookup[kv[0]] = kv[1];
			id++;
		}
		itr++;
	}
	std::map<std::string, std::string>::iterator indNotFound = phenotypeLookup.end();
	std::vector<Individual>::iterator indItr = individuals.begin();
	std::vector<Individual>::iterator indEnd = individuals.end();
	while (indItr != indEnd) {
		if (phenotypeLookup.find(indItr->indID) != indNotFound)
			indItr->status = atof(phenotypeLookup[indItr->indID].c_str());
		indItr++;
	}
}
*/
/*


void BinApplication::GenerateBinContentLookup(std::multimap<uint, uint>& binContents) {
	binContents.clear();
	std::vector<std::set<uint> > contributors;
	GetBinContributors(contributors);
	
	
	std::map<uint, uint>::iterator itr = binIndex.begin();
	std::map<uint, uint>::iterator end = binIndex.end();
	while (itr != end) {
		std::set<uint>& contribs = contributors[itr->second];
		std::set<uint>::iterator citr = contribs.begin();
		std::set<uint>::iterator cend = contribs.end();
		while (citr != cend) {
			binContents.insert(std::make_pair(itr->first, *(citr)));
			if ((*citr) > dataset.Size())
				std::cerr<<" Oversized Index ("<<dataset.Size()<<"):\t"<<regions[itr->first].id<<"\t"<<regions[itr->first].name<<"\t"<<*citr<<"\n";
			citr++;
		}
		itr++;
	}
}
*/

}
