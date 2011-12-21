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

#include <sqlite3.h>

// Use the boost filesystem library to work with OS-independent paths
#include <boost/filesystem.hpp>
//using boost::filesystem;

// Rewriting regionManager and groupManager from scratch
//#include "knowledge/regionmanagerdb.h"
//#include "knowledge/groupmanagerdb.h"

//#include "knowledge/genegenemodelarchive.h"
//#include "knowledge/snpdataset.h"
//#include "knowledge/def.h"
//#include "liftover/converterdb.h"

using std::string;
using std::vector;
using std::map;

#include "knowledge/Locus.h"
#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Information.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/GroupCollectionSQLite.h"

using std::string;

namespace BioBin {

class Application {
public:
	
	Application();
	virtual ~Application();

	void Init(const string& dbFilename, bool reportVersions);

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
	uint LoadGroupDataByName(T1_cont& userDefinedGroups, const set<uint>& group_type_ids);

	uint GetPopulationID(const string& pop);

	Knowledge::RegionCollection* GetRegions();

	/**
	 * Streams population IDs to stdout
	 */
	void ListPopulationIDs(std::ostream& os);

	/**
	 * Generates group/ID report showing meta groups, their children and all IDs in the current dataset
	 */
	void ListGroupIDs(std::ostream& os, const vector<string>& searchList);

	/**
	 * Update the database with new variation filename (allow users to move the file to a permanent location specific to their system)
	 */
	//void SetVariationFilename(const string& variationFilename);
	void ListGenes(std::ostream& os, const vector<string>& aliasList, const vector<string>& aliasType);

	void SetGeneExtension(uint geneBoundaryExt);
	void SetReportPrefix(const string& pref);
	void UseHtmlReports(bool useHtml);

	void LoadBuildConverter(const string& build);

	// Returns a vector of the category (source) IDs
	vector<uint> ManagerIDs();
	Knowledge::GroupCollection* GroupManager(uint idx);

	static bool errorExit;										///< When exiting on errors, we won't report the files that "would" have been generated.

protected:
	uint GetPopID(const string& pop);

	///< The name of the database file
	std::string dbFilename;

	///< The database session - going to not use this
	//soci::session sociDB;
	// Instead, use a standard sqlite3 database.
	sqlite3* _db;

	// Basic information about the database
	Knowledge::Information* _info;

	///< The genes
	Knowledge::RegionCollection* regions;
	///< The knowedge meta groups
	std::map<uint, Knowledge::GroupCollection*> groups;
	///< the data associated with the user
	vector<Knowledge::Locus*> dataset;

	///< The variation version (to guarantee that the variations file is correct for the database being used)
	uint varVersion;
	///< Filename for variation data-this might not be opened, depending on how the user loads their data
	std::string variationFilename;
	///< The list of report filenames generated
	std::stringstream reportLog;

	///< Length of extension on either side of a gene (not to be mixed with LD extension)
	uint geneExtensionLength;
	///< Turn on/off HTML report generation
	bool htmlReports;
	///< Report Prefix
	std::string reportPrefix;

	///< Converter structure for map files (using liftover chains)
	Knowledge::Liftover::Converter* buildConverter;

	//std::multimap<uint, uint> geneLookup;

private:




};


template <class T_cont>
uint Application::LoadRegionData(T_cont& aliasesNotFound, const vector<string>& aliasList) {

	regions->Load(aliasList);
	regions->associateLoci(dataset.begin(), dataset.end());

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
uint Application::LoadGroupDataByName(T1_cont& userDefinedGroups,
		const set<uint>& group_type_ids) {

	map<int, string> group_types;

	_info->getGroupTypes(group_type_ids, group_types);

	map<int, string>::const_iterator itr = group_types.begin();
	map<int, string>::const_iterator end = group_types.end();

	uint totalGroupsLoaded = 0;
	int max_id = 0;
	int curr_id;
	while(itr != end){
		Knowledge::GroupCollection* new_group = (Knowledge::GroupCollection*)
				new Knowledge::GroupCollectionSQLite((*itr).first, (*itr).second, _db);
		new_group->Load(*regions);
		curr_id = (*itr).first;
		groups[curr_id] = new_group;
		totalGroupsLoaded += new_group->size();
		max_id = curr_id > max_id ? curr_id : max_id;
		++itr;
	}

	vector<string> unmatchedAliases;
	typename T1_cont::const_iterator udItr = userDefinedGroups.begin();
	typename T1_cont::const_iterator udEnd = userDefinedGroups.end();

	while (udItr != udEnd) {
		string fn = *udItr;
		//Give some bogus groupType, since it will be found in the file
		Knowledge::GroupCollection* new_group =
				new Knowledge::GroupCollectionSQLite(++max_id, fn, _db);
		totalGroupsLoaded +=
				new_group->LoadArchive(*regions, fn, unmatchedAliases);
		groups[max_id] = new_group;
		++udItr;
	}

	return totalGroupsLoaded;

}


}


#endif	/* APPLICATION_H */

