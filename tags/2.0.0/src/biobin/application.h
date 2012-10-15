/* 
 * File:   Application.h
 * Author: torstees
 *
 * Created on March 28, 2011, 1:45 PM
 */

#ifndef APPLICATION_H
#define	APPLICATION_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <new>
#include <deque>

#include <sqlite3.h>

// Use the boost filesystem library to work with OS-independent paths
#include <boost/filesystem.hpp>

#include "knowledge/Locus.h"
#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Information.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/GroupCollectionSQLite.h"

using std::string;
using std::vector;
using std::map;
using std::new_handler;
using std::set_new_handler;
using std::deque;

namespace BioBin {

class Application {

public:
	
	Application(const string& db_fn);
	virtual ~Application();

	string GetReportLog();
	string AddReport(const string& suffix, const string& extension, const string& description);

	// We want to separate the next two steps because some jobs require only one
	// of the two be performed, and they do take time to complete
	/**
	 * Loads region data from database
    * @param pop indicates LD boundaries to be used, if any (if NO-LD, then the variable geneExtensionLength is used for static extensions
    * @param regionAliases A list of region aliases to be loaded. This will restrict the genes loaded to only those found (this is intended for certain types of reporting only)
    * @return Number of regions loaded
    */
	template <class T_cont>
	uint LoadRegionData(T_cont& aliasesNotFound,
			const vector<string>& aliasList);		///< Comma separated list of strings

	/**
	 * Loads group data using group name rather than IDs
    * @param userDefinedGroups
    * @param groupNames
    * @return
    */
	template <class T1_cont>
	uint LoadGroupDataByName(T1_cont& userDefinedGroups);

	uint GetPopulationID(const string& pop);

	Knowledge::RegionCollection* GetRegions();

	void SetGeneExtension(uint geneBoundaryExt);
	void SetReportPrefix(const string& pref);

	static bool errorExit;										///< When exiting on errors, we won't report the files that "would" have been generated.
	static std::string reportPrefix;
	static bool c_print_sources;
	static bool c_print_populations;
	static bool s_run_normal;

private:
	void Init(const string& dbFilename, bool reportVersions);

protected:
	uint GetPopID(const string& pop);

	///< The name of the database file
	std::string dbFilename;

	// Instead, use a standard sqlite3 database.
	sqlite3* _db;

	// Basic information about the database
	Knowledge::Information* _info;

	///< The genes
	Knowledge::RegionCollection* regions;
	///< The knowedge meta groups
	Knowledge::GroupCollection* groups;
	///< the data associated with the user
	deque<Knowledge::Locus*> dataset;

	///< The variation version (to guarantee that the variations file is correct for the database being used)
	uint varVersion;
	///< Filename for variation data-this might not be opened, depending on how the user loads their data
	std::string variationFilename;
	///< The list of report filenames generated
	std::stringstream reportLog;

	///< Length of extension on either side of a gene (not to be mixed with LD extension)
	uint geneExtensionLength;

//Everything from here on down has to do with installing a new handler that
// will try to get sqlite to give up some of its cache
public:
	// A function to release as much db cache as possible
	static void releaseDBCache();

private:

	static new_handler currentHandler;
};


template <class T_cont>
uint Application::LoadRegionData(T_cont& aliasesNotFound, const vector<string>& aliasList) {

	regions->Load(aliasList);

	//T_cont::const_iterator pos = aliasesNotFound.end();

	vector<string>::const_iterator itr = aliasList.begin();
	vector<string>::const_iterator end = aliasList.end();
	while (itr != end) {
		if (!regions->isValid(*itr))
			aliasesNotFound.insert(aliasesNotFound.end(), *itr);
		++itr;
	}

	return aliasList.size() - aliasesNotFound.size();
}


/**
 * This function depends on regions having been properly loaded
 */
template <class T1_cont>
uint Application::LoadGroupDataByName(T1_cont& userDefinedGroups) {

	groups = new Knowledge::GroupCollectionSQLite(*regions, _db, _info);
	groups->Load();

	vector<string> unmatchedAliases;
	typename T1_cont::const_iterator udItr = userDefinedGroups.begin();
	typename T1_cont::const_iterator udEnd = userDefinedGroups.end();

	while (udItr != udEnd) {
		groups->LoadArchive(*udItr, unmatchedAliases);
		++udItr;
	}


	return groups->size();

}


}


#endif	/* APPLICATION_H */

