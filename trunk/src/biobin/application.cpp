/* 
 * File:   Application.cpp
 * Author: torstees
 * 
 * Created on March 28, 2011, 1:45 PM
 */

#include "application.h"

#include <iomanip>
#include <sstream>
#include <iostream>

#include "knowledge/liftover/ConverterSQLite.h"
#include "knowledge/InformationSQLite.h"
#include "knowledge/RegionCollectionSQLite.h"
#include "knowledge/GroupCollectionSQLite.h"

namespace BioBin {

bool Application::errorExit = false;
std::string Application::reportPrefix = "biobin";
bool Application::c_print_sources = false;
bool Application::c_print_populations = false;
bool Application::s_run_normal = true;

new_handler Application::currentHandler;

Application::Application(const string& db_fn) :
		dbFilename(db_fn), varVersion(0), geneExtensionLength(0) {

	Init(db_fn, true);

	if(c_print_populations){
		_info->printPopulations(std::cout);
	}

	if(c_print_sources){
		_info->printSources(std::cout);
	}

	if(!currentHandler){
		currentHandler = set_new_handler(releaseDBCache);
	}
}

Application::~Application(){

	sqlite3_close(_db);

	if (!errorExit){
		std::cerr<<GetReportLog()<<"\n";
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

Knowledge::RegionCollection* Application::GetRegions() {
	return regions;
}

void Application::SetGeneExtension(uint geneBoundaryExt) {
	geneExtensionLength = geneBoundaryExt;
}

void Application::SetReportPrefix(const string& pref) {
	if (pref=="")
		reportPrefix = "biobin";
	else
		reportPrefix = pref;
}

string Application::GetReportLog() {
	return reportLog.str();
}

string Application::AddReport(const string& suffix, const string& extension, const string& description) {
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

uint Application::GetPopulationID(const string& pop_str) {
	return _info->getPopulationID(pop_str);
}

void Application::Init(const string& filename, bool reportVersion) {
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
	_info = new Knowledge::InformationSQLite(_db);
	regions = new Knowledge::RegionCollectionSQLite(_db, dataset, _info);

}

void Application::releaseDBCache(){
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
