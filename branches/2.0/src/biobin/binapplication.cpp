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

BinApplication::BinApplication(const string& vcf_file) :
		binData(_pop_mgr), _data(vcf_file) {}

void BinApplication::InitBins() {
	
	binData.InitBins(groups, *regions, dataset, _info);

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
			++itr;
		}
	}
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

void BinApplication::writeLoci(const string& filename, const string& sep) const{
	std::ofstream locusFile(filename.c_str());
	locusFile << "Chromosome" << sep << "Location" << sep << "ID" << sep
			<< "Alleles" << sep << "Case Allele Freq." << sep << "Rare" << sep
			<< "gene(s)" << sep << "bin name(s)\n";
	vector<Knowledge::Locus*>::const_iterator itr = dataset.begin();
	vector<Knowledge::Locus*>::const_iterator end = dataset.end();

	unordered_map<Knowledge::Locus*, float> case_maf;

	_data.getCaseAF(dataset, _pop_mgr.getControls(), case_maf);

	while(itr != end){
		(*itr)->print(locusFile, sep);
		locusFile << sep;
		(*itr)->printAlleles(locusFile, "|");
		locusFile << sep << case_maf[*itr] << sep
				<< static_cast<int>((*itr)->isRare()) << sep;
		// Print the genes here
		Knowledge::RegionCollection::const_region_iterator r_itr =
				regions->positionBegin((*itr)->getChrom(), (*itr)->getPos());
		Knowledge::RegionCollection::const_region_iterator r_end =
				regions->positionEnd((*itr)->getChrom(), (*itr)->getPos());
		if (r_itr != r_end){
			locusFile << (*r_itr)->getName();
			while(++r_itr != r_end){
				locusFile << "|" << (*r_itr)->getID();
			}
		}
		locusFile << sep;
		// Now print the bins
		binData.printBins(locusFile, *itr, "|");
		locusFile << "\n";
		++itr;
	}
	locusFile.close();
}

void BinApplication::writeAFData(const string& filename, const string& sep) const{
	std::ofstream freqFile(filename.c_str());

	unordered_map<Knowledge::Locus*, float> case_maf;

	_data.getCaseAF(dataset, _pop_mgr.getControls(), case_maf);

	vector<Knowledge::Locus*>::const_iterator itr = dataset.begin();
	vector<Knowledge::Locus*>::const_iterator end = dataset.end();

	freqFile << "Locus" << sep << "Control NMAF" << sep << "Case NMAF" << sep
			<< "Rare" << sep << "Bins\n";
	while(itr != end){
		freqFile << (*itr)->getID() << sep << 1 - (*itr)->majorAlleleFreq() << sep
				<< case_maf[*itr] << sep
				<< static_cast<int>((*itr)->isRare()) << sep;

		binData.printBins(freqFile, *itr, "|");

		freqFile << "\n";

		++itr;
	}

	freqFile.close();

}

void BinApplication::writeBinFreqData(const string& filename, const string& sep) const {
	std::ofstream freqFile(filename.c_str());
	_pop_mgr.printBinFreq(freqFile, binData, sep);
	freqFile.close();

}

}
