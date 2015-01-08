/* 
 * File:   binapplication.cpp
 * Author: torstees
 * 
 * Created on June 22, 2011, 10:35 AM
 */

#include "binapplication.h"

#include <iomanip>

#include <boost/algorithm/string.hpp>

#include "knowledge/InformationSQLite.h"
#include "knowledge/RegionCollectionSQLite.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::deque;
using std::new_handler;
using std::set_new_handler;
using std::stringstream;
using std::ostream;

using boost::unordered_map;

namespace BioBin {

bool BinApplication::c_transpose_bins = false;
bool BinApplication::errorExit = false;
std::string BinApplication::reportPrefix = "biobin";
bool BinApplication::c_print_sources = false;
bool BinApplication::c_print_populations = false;
bool BinApplication::s_run_normal = true;

new_handler BinApplication::currentHandler;

BinApplication::BinApplication(const string& db_fn, const string& vcf_file) :
	dbFilename(db_fn), varVersion(0), geneExtensionLength(0),
			_pop_mgr(vcf_file), binData(_pop_mgr) {

	Init(db_fn, true);

	if (c_print_populations) {
		_info->printPopulations(std::cout);
	}

	if (c_print_sources) {
		_info->printSources(std::cout);
	}

	if (!currentHandler) {
		currentHandler = set_new_handler(releaseDBCache);
	}
}

BinApplication::~BinApplication(){

	sqlite3_close(_db);

	if (!errorExit){
		std::cout<<GetReportLog()<<"\n";
	}

	delete _info;
	delete regions;

	deque<Knowledge::Locus*>::iterator d_itr = dataset.begin();
	while(d_itr != dataset.end()){
		delete *d_itr;
		++d_itr;
	}
	//dataset.clear();

	delete groups;
}

void BinApplication::InitBins() {
	
	_info->loadWeights(*regions);

	binData.InitBins(*regions, dataset, _info);

	std::cout<<"\n   Total SNPS:   "<<std::setw(10)<<std::right<<dataset.size()<<"\n"
				<<"   Variants:     "<<std::setw(10)<<std::right<<binData.numVariants()<<"\n"
				<<" * Rare Variants:"<<std::setw(10)<<std::right<<binData.numRareVariants()<<"\n"
				<<"   Total Bins:   "<<std::setw(10)<<std::right<<binData.numBins()<<"\n";

	std::cout<<"\n   * Rare variants are those whose minor allele frequency is below "
			 <<BinManager::mafCutoff<<" and above "<<BinManager::mafThreshold<<"\n\n";

	if (binData.numBins() > 0 && binData.numBins() < 500) {
		std::cout<<"\n\nBin Name\tVariant Count\n";
		BinManager::const_iterator itr = binData.begin();
		BinManager::const_iterator end = binData.end();
		while(itr != end){
			std::cout << (*itr)->getName() << "\t" << (*itr)->getSize() << "\n";
			++itr;
		}
	}

	std::cout << std::endl;
	binData.printLocusBinCount(std::cout);
	std::cout << std::endl;

	// Add output about number of bins
}

void BinApplication::writeBinData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	if(c_transpose_bins){
		_pop_mgr.printBinsTranspose(file, binData, *_info, sep);
	}else{
		_pop_mgr.printBins(file, binData, *_info, sep);
	}
	file.close();
}

/*void BinApplication::writeGenotypeData(const string& filename, const string& sep) const{
	std::ofstream file(filename.c_str());
	_pop_mgr.printGenotypes(file, dataset, sep);
	file.close();
}
*/

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

/*void BinApplication::writeAFData(const string& filename, const string& sep) const{
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

		freqFile << sep << _pop_mgr.getControlAF(**itr) << sep
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
*/

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

string BinApplication::AddReport(const string& suffix, const string& extension, const string& description) {
	string newFilename(reportPrefix);
	if (suffix.size() > 0){
		newFilename += ("-" + suffix);
	}
	if (extension.size() > 0){
		newFilename += ("." + extension);
	}
	reportLog<<std::setw(50)<<std::right<<newFilename<<" : "<<description<<"\n";
	return newFilename;
}

void BinApplication::Init(const string& filename, bool reportVersion) {
	dbFilename = filename;
	boost::filesystem::path dbPath = boost::filesystem::path(dbFilename);
	bool fileFound = false;
	if (boost::filesystem::is_regular_file(dbPath)) {
		fileFound = true;
	}else{
		#ifdef DATA_DIR
			if (dbPath.is_relative()){
				dbPath = (boost::filesystem::path(std::string(DATA_DIR))/=(dbPath));
				if (boost::filesystem::is_regular_file(dbPath)){
					fileFound=true;
				}
			}
		#endif
	}

	if (!fileFound){
		throw std::ios_base::failure("File " + filename + " not found");
	}

	sqlite3_open(dbPath.c_str(), &_db);

	// set the pragma to only allow temporary storage in memory
	string memory_pragma = "PRAGMA temp_store=2;";
	sqlite3_exec(_db, memory_pragma.c_str(), NULL, NULL, NULL);

	_info = new Knowledge::InformationSQLite(_db);
	regions = new Knowledge::RegionCollectionSQLite(_db, dataset, _info);

}

void BinApplication::releaseDBCache(){
	// If we are here, we are damn near out of memory, so try to get SQLite
	// to release some memory (10MB if possible)!
	int mem_released = sqlite3_release_memory(10000);

	// If we didn't release any memory, abandon hope!
	if(!mem_released){
		set_new_handler(currentHandler);
		currentHandler = 0;
	}
}

}
