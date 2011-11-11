/* 
 * File:   binapplication.cpp
 * Author: torstees
 * 
 * Created on June 22, 2011, 10:35 AM
 */

#include "binapplication.h"

namespace BioBin {

BinApplication::BinApplication() {
}


BinApplication::~BinApplication() {
}


void BinApplication::InitVcfDataset(std::string& filename, std::string& genomicBuild, Knowledge::SnpDataset& lostSnps, std::vector<uint>& locusRemap, DataImporter& vcfimporter) {
	
	std::cerr<<"Loading VCF Data\n";
	if (vcfimporter.Open(filename.c_str(), (char)-1)) {
		std::vector<Utility::Locus> locusArray;
		vcfimporter.GetAllAlleleFrequencies(locusArray);
		locusRemap.clear();
		locusRemap.reserve(locusArray.size());
		//const std::vector<Utility::Locus>& locusArray = vcfimporter.GetLoci();
		LiftOver::ConverterDB cnv;
		int chainCount = cnv.LoadFromDB(genomicBuild.c_str(), sociDB);
		
		if (chainCount > 0) {
			std::string conversionLog = this->AddReport("lift-over", "tsv", "SNPs that were lifted over to new build which differed dramatically or changed chromosome");
			std::ofstream cnvLog(conversionLog.c_str());
			cnvLog<<"RSID,Chrom(Orig),Pos(Orig),Chrom(New),Pos(New)\n";
			
			std::multimap<Utility::Locus, Utility::Locus> converted;
			std::multimap<Utility::Locus, Utility::Locus>::iterator missing;
			cnv.ConvertDataset(locusArray, converted);
			Utility::SnpArray::iterator itr = locusArray.begin();
			Utility::SnpArray::iterator end = locusArray.end();

			uint i=0;
			//uint validLocus = 0;
			std::stringstream missingSNPs;
			while (itr != end) {
				Utility::Locus &orig = *itr;
				
				if (converted.find(orig) != missing) {
					std::multimap<Utility::Locus, Utility::Locus>::iterator first = converted.lower_bound(orig);
					std::multimap<Utility::Locus, Utility::Locus>::iterator last  = converted.upper_bound(orig);
					
					if (converted.count(orig) != 1)
						std::cerr<<"It was observed that there are multiple hits returned by convert dataset: "<<orig.RSID()<<" has "<<converted.count(orig)<<" counterparts.\n";
					while (first != last) {
						if (first->second.pos == 0)
							missingSNPs<<first->first.RSID()<<"\t"<<Utility::ChromFromIntChr(first->first.chrom)<<"\t"<<first->first.pos<<"\n";
						else {
							if (first->second.chrom > 0) {
								if (first->first.Chrom() != first->second.Chrom() || ((float)abs((float)first->first.pos - (float)first->second.pos)/(float)first->first.pos)> 0.01)
									cnvLog<<first->first.RSID()<<"\t"
										<<first->first.Chrom()<<"\t"
										 <<first->first.pos<<"\t"
										 <<first->second.Chrom()<<"\t"
										 <<first->second.pos<<"\t"
										 <<first->second.RSID()<<"\n";
								dataset.AddSNP(first->second);
								//dataset.AddSNP(first->second.chrom, first->second.pos, first->second.RSID().c_str());
								//locusRemap[first->second.chrom].push_back(validLocus++);
								locusRemap.push_back(i);
								locusArray[i] = first->second;
							}		
							else 
								cnvLog<<first->first.RSID()<<"\t"
									<<first->first.Chrom()<<"\t"
									 <<first->first.pos<<"\t"
									 <<" **** \n";
						}
						first++;
					}
				}
				itr++;
				i++;
			}
			if (missingSNPs.str().length() > 0) {
				std::string filename = AddReport("missing-snps", "txt", "SNPs that were dropped during build conversion");
				std::ofstream file(filename.c_str());
				file<<missingSNPs.str();
			}
		} else {
			uint locusCount = locusArray.size();
			dataset.AddSNPs(locusArray);
			for (uint i=0; i<locusCount; i++) 
				locusRemap.push_back(i);	
		}
	}
}
//std::vector<uint> locusRemap;						///< Necessary to translate from local indexes to those in the vcf file


std::pair<uint, uint> BinApplication::InitBins(std::vector<uint>& locusRemap, 
			DataImporter& vcfimporter) {
	
	Utility::IdCollection variants;
	Utility::IdCollection rareVariants;

	//if (vcfimporter.Open(filename.c_str(), (char)-1)) {
		//We now have our snp dataset set up-so it's time to start the binning process
	std::pair<uint, uint> binGenotypeCounts;
	binGenotypeCounts = binData.InitBins(groups, regions, dataset, locusRemap);
	binData.CollectVariantGroups(variants, rareVariants);
	std::cerr<<"   Total SNPS:   "<<std::setw(10)<<std::right<<dataset.Size()<<"\n"
				<<"   Variants:     "<<std::setw(10)<<std::right<<variants.size()<<"\n"
				<<" * Rare Variants:"<<std::setw(10)<<std::right<<rareVariants.size()<<"\n"
				<<"   Total Bins:   "<<std::setw(10)<<std::right<<binGenotypeCounts.first<<"\n";

	std::cerr<<"\n   * Rare variants are those whose minor alleles sum is below: "<<BinManager::mafCutoff<<"\n";

	if (binGenotypeCounts.first < 500) {
		std::cerr<<"\n\nBin Name\tSNP Count\n";
		std::vector<std::set<uint> > contributors;
		binData.BuildContributorList(contributors);
		uint count = contributors.size();
		for (uint i=0; i<count; i++) 
			std::cerr<<binData.BinName(i)<<"\t"<<contributors[i].size()<<"\n";

	}
	
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

	uint locusCount = dataset.Size();
	std::string ofn = AddReport("locus", "csv", "Locus Description");
	std::ofstream locusFile(ofn.c_str());
	locusFile<<"Chromosome,Location,ID,all(1):freq(1),all(2):freq(2),type,gene(s),bin name(s)\n";
	for (uint i=0; i<locusCount; i++) {
		dataset[i].Print(locusFile, ",");
		binData.DescribeLocus(i, locusFile, regions,dataset);
	}
	locusFile.close();

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
}


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
}
