/* 
 * File:   Application.cpp
 * Author: torstees
 * 
 * Created on March 28, 2011, 1:45 PM
 */

#include "application.h"

#include <iomanip>
#include <sstream>

#include <sqlite3.h>

#include "knowledge/liftover/ConverterSQLite.h"
#include "knowledge/InformationSQLite.h"
#include "knowledge/RegionCollectionSQLite.h"
#include "knowledge/GroupCollectionSQLite.h"

namespace BioBin {

bool Application::errorExit = false;
std::string Application::reportPrefix = "biobin";

Application::Application() : dbFilename(""), varVersion(0), geneExtensionLength(0), htmlReports(false) {}

Application::~Application(){
	sqlite3_close(_db);

	if (!errorExit){
		std::cerr<<GetReportLog()<<"\n";
	}

	delete _info;
	delete regions;

	vector<Knowledge::Locus*>::iterator d_itr = dataset.begin();
	while(d_itr != dataset.end()){
		delete *d_itr;
		++d_itr;
	}
	dataset.clear();

	map<uint, Knowledge::GroupCollection*>::iterator g_itr = groups.begin();
	while(g_itr != groups.end()){
		delete (*g_itr).second;
		++g_itr;
	}
	groups.clear();


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

void Application::UseHtmlReports(bool doUse) {
	htmlReports = doUse;
}

vector<uint> Application::ManagerIDs() {
	vector<uint> ids;

	map<uint, Knowledge::GroupCollection*>::const_iterator itr = groups.begin();
	map<uint, Knowledge::GroupCollection*>::const_iterator end = groups.end();
	while (itr != end)
		ids.push_back(itr++->first);
	return ids;
}

Knowledge::GroupCollection* Application::GroupManager(uint idx) {
	map<uint, Knowledge::GroupCollection*>::const_iterator itr = groups.find(idx);
	map<uint, Knowledge::GroupCollection*>::const_iterator end = groups.end();

	// Return the pointer to GroupManager if found, o/w return null pointer
	return ((!(itr == end)) ? (*itr).second : (NULL));
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
	regions = new Knowledge::RegionCollectionSQLite(_db, dataset);

}



}
