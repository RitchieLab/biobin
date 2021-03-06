/* 
 * File:   binapplication.cpp
 * Author: torstees
 * 
 * Created on June 22, 2011, 10:35 AM
 */

#include "binapplication.h"

#include "main.h"

#include <iomanip>
#include <cstdio>

#include <boost/algorithm/string.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>

#include "knowledge/InformationSQLite.h"
#include "knowledge/RegionCollectionSQLite.h"

#include "util/Phenotype.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::deque;
using std::new_handler;
using std::set_new_handler;
using std::stringstream;
using std::ostream;
using std::istream;

using boost::unordered_map;

namespace BioBin {

std::string BinApplication::reportPrefix = "biobin";
bool BinApplication::c_transpose_bins = false;
bool BinApplication::errorExit = false;
bool BinApplication::c_print_sources = false;
bool BinApplication::c_print_populations = false;
bool BinApplication::s_run_normal = true;
unsigned int BinApplication::n_threads = 0;

new_handler BinApplication::currentHandler;

BinApplication::BinApplication(const string& db_fn, const string& vcf_file) :
	dbFilename(db_fn), _info(0), regions(0), groups(0), varVersion(0), geneExtensionLength(0),
			_pop_mgr(vcf_file) {

	Init(db_fn, true);

	_pop_mgr.setInfo(_info);

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

	if(_info){
		delete _info;
	}

	if(groups){
		delete groups;
	}

	if(regions){
		delete regions;
	}

	deque<Knowledge::Locus*>::iterator d_itr = dataset.begin();
	while(d_itr != dataset.end()){
		delete *d_itr;
		++d_itr;
	}
	//dataset.clear();

	for(unsigned int i=0; i<_locus_bins.size(); i++){
		fclose(_locus_bins[i]);
	}
	_locus_bins.clear();

}

void BinApplication::binPhenotypes(PopulationManager::const_pheno_iterator& ph_itr){

	_pheno_mutex.lock();
		while(ph_itr != _pop_mgr.endPheno()){
			// use the default copy constructor for Utility::Phenotype
			Utility::Phenotype ph(*ph_itr);
			++ph_itr;
			_pheno_mutex.unlock();

			BinManager binData(_pop_mgr, *regions, dataset, *_info, ph);

			_output_mutex.lock();
			std::cout << "Phenotype: " << _pop_mgr.getPhenotypeName(ph.getIndex()) << std::endl;

			std::cout<<"\n   Loci:     "<<std::setw(10)<<std::right<<dataset.size()<<"\n"
						<<" * Rare Loci:"<<std::setw(10)<<std::right<<binData.numRareVariants()<<"\n"
						<<"   Total Bins:   "<<std::setw(10)<<std::right<<binData.size()<<"\n";

			std::cout<<"\n   * Rare variants are those whose minor allele frequency is below "
					 <<BinManager::mafCutoff<<" and above "<<BinManager::mafThreshold<<"\n\n";

			if (binData.size() > 0 && binData.size() < 100) {
				std::cout<<"\n\nBin Name\tVariant Count\tLoci Count" << std::endl;
				BinManager::const_iterator itr = binData.begin();
				BinManager::const_iterator end = binData.end();
				while(itr != end){
					std::cout << (*itr)->getName() << "\t" << (*itr)->getSize() << "\t" << (*itr)->getVariantSize() << "\n";
					++itr;
				}
			}

			std::cout << std::endl;
			binData.printLocusBinCount(std::cout);
			std::cout << std::endl;

			_output_mutex.unlock();

			// If we are creating a locus report, print all bin data for each locus
			if(Main::WriteLociData){
				FILE* fp = binData.printLocusBins(dataset,
						_pop_mgr.getPhenotypeName(ph.getIndex()));

				_data_mutex.lock();
				_locus_bins[ph.getIndex()] = fp;
				_data_mutex.unlock();
			}

			// print the Bin data
			if(Main::WriteBinData){
				std::string phenoname = _pop_mgr.getPhenotypeName(ph.getIndex());
				if(phenoname.size() > 0){
					phenoname = phenoname + "-";
				}

				std::string filename = reportPrefix + "-" + phenoname + "bins.csv";
				std::ofstream file(filename.c_str());
				binData.printBinData(file, Main::OutputDelimiter, c_transpose_bins);
				file.close();
			}

			_pheno_mutex.lock();
		}
		_pheno_mutex.unlock();

}

void BinApplication::InitBins() {

	_info->loadWeights(*regions);

	if(Main::WriteLociData){
		_locus_bins.resize(_pop_mgr.getNumPhenotypes(), 0);
	}
	PopulationManager::const_pheno_iterator ph_itr = _pop_mgr.beginPheno();

	//TODO: make this a threaded call
	if(n_threads == 0){
		binPhenotypes(ph_itr);
	} else {
		boost::thread_group tg;
		for (unsigned int i = 0; i < n_threads; i++) {
			tg.create_thread(boost::bind(&BinApplication::binPhenotypes, this, boost::ref(ph_itr)));
		}
		tg.join_all();
	}
}

void BinApplication::InitVcfDataset(const std::string& genomicBuild) {

	Knowledge::Liftover::ConverterSQLite cnv(genomicBuild, _db);
	int chainCount = cnv.Load();

	if (chainCount > 0) {
		_pop_mgr.loadLoci(dataset,reportPrefix,Main::OutputDelimiter,&cnv);
	}else{
		_pop_mgr.loadLoci(dataset,reportPrefix,Main::OutputDelimiter);
	}

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
	printEscapedString(locusFile, "Gene(s)", sep, sep_repl);

	vector<istream*> locus_streams;
	locus_streams.reserve(_locus_bins.size());
	string pheno_bins;
	for(unsigned int i=0; i<_locus_bins.size(); i++){
		rewind(_locus_bins[i]);
		locus_streams.push_back(
				new boost::iostreams::stream<boost::iostreams::file_descriptor>(
						boost::iostreams::file_descriptor(fileno(_locus_bins[i])),
						std::ios_base::binary | std::ios_base::in)
				);
		locusFile << sep;
		getline(*locus_streams[i], pheno_bins);
		printEscapedString(locusFile, pheno_bins + " Bin Name(s)", sep, sep_repl);
	}
	locusFile << "\n";

	deque<Knowledge::Locus*>::const_iterator itr = dataset.begin();


	while(itr != dataset.end()){
		printEscapedString(locusFile, (*itr)->getChromStr(), sep, sep_repl);
		locusFile << sep << (*itr)->getPos() << sep;
		printEscapedString(locusFile, (*itr)->getID(), sep, sep_repl);
		locusFile << sep;

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

		// print all of the old bins now
		for(unsigned int i=0; i<locus_streams.size(); i++){
			locusFile << sep;
			getline(*locus_streams[i], pheno_bins);
			printEscapedString(locusFile, pheno_bins, sep, sep_repl);
		}

		locusFile << "\n";
		++itr;
	}
	locusFile.close();
	for(unsigned int i=0; i<locus_streams.size(); i++){
		delete locus_streams[i];
	}
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
