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
//#include <soci-sqlite3.h>
//#include "utility/filetools.h"
//#include "taskgenecoverage.h"
//#include "ldsplineimporter.h"

#include "knowledge/liftover/ConverterSQLite.h"
#include "knowledge/InformationSQLite.h"
#include "knowledge/RegionCollectionSQLite.h"
#include "knowledge/GroupCollectionSQLite.h"
//#include "liftover/converterdb.h"

#include "utility/exception.h"

namespace BioBin {

bool Application::errorExit = false;

Application::Application() : dbFilename(""), varVersion(0), geneExtensionLength(0), htmlReports(false), reportPrefix("report") {}

Application::~Application(){
	sqlite3_close(_db);

	if (!errorExit){
		std::cerr<<GetReportLog()<<"\n";
	}

	delete _info;
	delete regions;
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

void Application::ListGenes(std::ostream& os, const vector<string>& aliasList, const vector<string>& aliasType) {
	_info->listRegions(os, aliasList, aliasType);
}

void Application::ListGroupIDs(std::ostream& os, const vector<string>& searchList) {
	_info->listGroupIDs(os, searchList);
}

void Application::ListPopulationIDs(std::ostream& os) {
	_info->listPopulationIDs(os);
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

void Application::LoadBuildConverter(const string& build) {
	buildConverter = new Knowledge::Liftover::ConverterSQLite(build, _db);
	buildConverter->Load();
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
		throw Utility::Exception::FileNotFound(filename.c_str());
	}

	sqlite3_open(filename.c_str(), &_db);
	_info = new Knowledge::InformationSQLite(_db);
	regions = new Knowledge::RegionCollectionSQLite(_db);

	string dbSnp = _info->getResourceVersion("ncbi");
	string ensembl = _info->getResourceVersion("ensembl");
	string hapmap = _info->getResourceVersion("hapmap");
	variationFilename = _info->getResourceVersion("variations");
	string build = _info->getResourceVersion("build");
	string variations;

	// TODO: I think this will always give 0 because we haven't even set the
	// variations string yet!!
	this->varVersion = atoi(variations.c_str());

	if (reportVersion) {
		std::cerr<<"\n------------------------- Dependency Versions ----------\n";
		std::cerr<<std::setw(38)<<std::right<<"dbSNP : "<<dbSnp<<"\n";
		std::cerr<<std::setw(38)<<std::right<<"Ensembl : "<<ensembl<<"\n";
		std::cerr<<std::setw(38)<<std::right<<"Hap Map LD : "<<hapmap<<"\n";
		std::cerr<<std::setw(38)<<std::right<<"Variation Filename : "<<variationFilename<<"\n";
		std::cerr<<std::setw(38)<<std::right<<"Genome Build : "<<build<<"\n";

	}

	// TODO: Load SNP role information here
	/*
	try {
		soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT id, role FROM snp_role");

		for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			soci::row const& row = *itr;
			uint id = row.get<int>(0);
			std::string name = row.get<std::string>(1);

			dataset.RoleDescription(id, name.c_str());
		}
	} catch (soci::soci_error const &e) {
		std::cerr<<"An error was encountered trying to load the role data. SNP role information will not be available.\n";
	}
	*/
}



}
