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

#include <boost/algorithm/string.hpp>

using std::string;
using std::vector;
using std::map;
using std::set;

using boost::unordered_map;

namespace BioBin {

bool BinApplication::c_transpose_bins = false;

BinApplication::BinApplication(const string& db_fn, const string& vcf_file) :
		Application(db_fn), _pop_mgr(vcf_file),	binData(_pop_mgr){}

void BinApplication::InitBins() {
	
	binData.InitBins(*regions, dataset, _info);

	std::cerr<<"\n   Total SNPS:   "<<std::setw(10)<<std::right<<dataset.size()<<"\n"
				<<"   Variants:     "<<std::setw(10)<<std::right<<binData.numVariants()<<"\n"
				<<" * Rare Variants:"<<std::setw(10)<<std::right<<binData.numRareVariants()<<"\n"
				<<"   Total Bins:   "<<std::setw(10)<<std::right<<binData.numBins()<<"\n";

	std::cerr<<"\n   * Rare variants are those whose minor alleles sum is below: "<<BinManager::mafCutoff<<"\n\n";

	if (binData.numBins() > 0 && binData.numBins() < 500) {
		std::cerr<<"\n\nBin Name\tVariant Count\n";
		BinManager::const_iterator itr = binData.begin();
		BinManager::const_iterator end = binData.end();
		while(itr != end){
			std::cerr << (*itr)->getName() << "\t" << (*itr)->getSize() << "\n";
			++itr;
		}
	}

	std::cerr << std::endl;
	binData.printLocusBinCount(std::cerr);
	std::cerr << std::endl;

	// Add output about number of bins
}

void BinApplication::writeBinData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	if(c_transpose_bins){
		_pop_mgr.printBinsTranspose(file, binData, sep);
	}else{
		_pop_mgr.printBins(file, binData, sep);
	}
	file.close();
}

void BinApplication::writeGenotypeData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	_pop_mgr.printGenotypes(file, dataset, sep);
	file.close();
}

void BinApplication::writeLoci(const string& filename, const string& sep) const{

	string sep_repl = getEscapeString(sep);

	std::ofstream locusFile(filename.c_str());

	printEscapedString(locusFile, "Chromosome", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Location", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "ID", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Alleles", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Case Allele Freq.", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Rare", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Gene(s)", sep, sep_repl);
	locusFile << sep;
	printEscapedString(locusFile, "Bin Name(s)", sep, sep_repl);
	locusFile << "\n";

	deque<Knowledge::Locus*>::const_iterator itr = dataset.begin();

	unordered_map<Knowledge::Locus*, float> case_maf;

	//_data.getCaseAF(dataset, _pop_mgr.getControls(), case_maf);

	while(itr != dataset.end()){
		printEscapedString(locusFile, (*itr)->getChromStr(), sep, sep_repl);
		locusFile << sep << (*itr)->getPos() << sep;
		printEscapedString(locusFile, (*itr)->getID(), sep, sep_repl);
		locusFile << sep;

		stringstream allele_str;
		(*itr)->printAlleles(allele_str, "|");
		printEscapedString(locusFile, allele_str.str(), sep, sep_repl);

		locusFile << sep << _pop_mgr.getCaseAF(**itr) << sep
				<< static_cast<int>((*itr)->isRare()) << sep;
		// Print the genes here
		Knowledge::RegionCollection::const_region_iterator r_itr =
				regions->locusBegin(*itr);
		Knowledge::RegionCollection::const_region_iterator r_end =
				regions->locusEnd(*itr);

		stringstream gene_str;
		if (r_itr != r_end){
			gene_str << (*r_itr)->getName();
			while(++r_itr != r_end){
				gene_str << "|" << (*r_itr)->getName();
			}
		}
		printEscapedString(locusFile, gene_str.str(), sep, sep_repl);
		locusFile << sep;
		// Now print the bins
		stringstream bin_str;
		binData.printBins(bin_str, *itr, "|");
		printEscapedString(locusFile, bin_str.str(), sep, sep_repl);
		locusFile << "\n";
		++itr;
	}
	locusFile.close();
}

void BinApplication::writeAFData(const string& filename, const string& sep) const{
	string sep_repl = getEscapeString(sep);

	std::ofstream freqFile(filename.c_str());

	unordered_map<Knowledge::Locus*, float> case_maf;

	//_data.getCaseAF(dataset, _pop_mgr.getControls(), case_maf);

	deque<Knowledge::Locus*>::const_iterator itr = dataset.begin();

	printEscapedString(freqFile, "Locus", sep, sep_repl);
	freqFile << sep;
	printEscapedString(freqFile, "Control NMAF", sep, sep_repl);
	freqFile << sep;
	printEscapedString(freqFile, "Case NMAF", sep, sep_repl);
	freqFile << sep;
	printEscapedString(freqFile, "Rare", sep, sep_repl);
	freqFile << sep;
	printEscapedString(freqFile, "Bins", sep, sep_repl);
	freqFile << "\n";

	while(itr != dataset.end()){
		printEscapedString(freqFile, (*itr)->getID(), sep, sep_repl);

		freqFile << sep << 1 - (*itr)->majorAlleleFreq() << sep
				<< _pop_mgr.getCaseAF(**itr) << sep
				<< static_cast<int>((*itr)->isRare()) << sep;

		stringstream ss;
		binData.printBins(ss, *itr, "|");
		printEscapedString(freqFile, ss.str(), sep, sep_repl);

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

void BinApplication::printEscapedString(ostream& os, const string& toPrint, const string& toRepl, const string& replStr) const{
	os << boost::algorithm::replace_all_copy(toPrint, toRepl, replStr);
}

string BinApplication::getEscapeString(const string& sep) const{
	string sep_repl = "_";
	if (sep == sep_repl){
		sep_repl = "-";
	}
	return sep_repl;
}

}
