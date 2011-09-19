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
#include "knowledge/regionmanagerdb.h"
#include "knowledge/groupmanagerdb.h"
#include "knowledge/genegenemodelarchive.h"
#include "knowledge/snpdataset.h"
#include "knowledge/def.h"
#include "knowledge/genegenemodelarchive.h"
#include "liftover/converterdb.h"

namespace Biofilter {

class Application {
public:
	
	Application();
	virtual ~Application();

	void InitBiofilter(const char *dbFilename, bool reportVersions);

	std::string GetReportLog();

	std::string AddReport(const char *suffix, const char *extension, const char *description);
	/**
	 * Loads data based on plink map file
    * @param filename The file used to specify the SNP data
	 * @param genomicBuild The build id associated with the map
	 * @param lostSnps array of SNPs who cound not be merged for some reason (because they were in regions that didn't properly lift over)
    * @return number of SNPs successfully loaded
    */
	uint LoadMapData(const char *filename, const char *genomicBuild, Knowledge::SnpDataset& lostSnps, bool doAlignData = false);

	/**
	 * Loads SNPs based on RS ID
    * @param filename The name of the file containing the RS IDs
    * @param lostSnps The SNPs that couldn't be found in the variations file
    * @return Number of SNPs successfully loaded
    */
	uint LoadSnpsSource(const char *filename, std::set<std::string>& lostSnps);

	// We want to separate the next two steps because some jobs require only one
	// of the two be performed, and they do take time to complete
	/**
	 * Loads region data from database
    * @param pop indicates LD boundaries to be used, if any (if NO-LD, then the variable geneExtensionLength is used for static extensions
    * @param regionAliases A list of region aliases to be loaded. This will restrict the genes loaded to only those found (this is intended for certain types of reporting only)
    * @return Number of regions loaded
    */
	uint LoadRegionData(const char *pop, Utility::StringArray& aliasesNotFound, Utility::StringArray& aliasList);		///< Comma separated list of strings
	/**
	 * Loads all group data (or subselection based on groupInclusions) as well as all user defined stuff
	 * @param userDefinedGroups Array of configuration strings to define the different groups
	 * @param groupInclusions The IDs associated for groups to be loaded (can include meta group ID)
    * @return Number of groups loaded
    */
	uint LoadGroupData(Utility::StringArray& userDefinedGroups, Utility::StringArray& groupInclusions);

	/**
	 * Loads group data using group name rather than IDs
    * @param userDefinedGroups
    * @param groupNames
    * @return
    */
	uint LoadGroupDataByName(Utility::StringArray& userDefinedGroups, Utility::StringArray& groupNames, Utility::IdCollection& ids);

	void GeneCoverage(Utility::StringArray& rsSources, Utility::StringArray& mapSources, const char *geneFilename, const char *population);

	uint GetPopulationID(const char *pop);

	void LoadLdSpline(const char *cfg);

	Knowledge::RegionManager *GetRegions();
	Knowledge::SnpDataset *GetDataset();
	std::multimap<uint, uint> *GetGeneLookup();
	Knowledge::GeneGeneModelArchive *GetGeneGeneModels();
	/**
	 * This can be used to allow users to provide group names, instead of group IDs 
	 * The return value should be sufficient to pass to LoadGroupData
	 */
	std::string ConvertGroupNamesToIDs(const char *groupList);

	/**
	 * Streams population IDs to stdout
	 */
	void ListPopulationIDs(std::ostream& os);

	/**
	 * Generates group/ID report showing meta groups, their children and all IDs in the current dataset
	 */
	void ListGroupIDs(std::ostream& os, Utility::StringArray& searchList);


	/**
	 * Generate simple alias type ID report
    */
	void ListAliasTypes(std::ostream& os);

	/**
	 * Produce gene reporting
    * If detailed, then SNPs will be included in the report
    */
	void ListGenes(std::ostream, bool detailed);

	std::multimap<uint, uint> BuildSnpGeneMap();

	/**
	 * Cheap way to turn add/remove indexes so that certain actions will be as fast as possible
	 */
	void StripOptimization();
	void PerformOptimization();

	/**
	 * Loads aliases based on the preferred type.
	 * @param aliasTypeList Comma separated list of type_ids (from the database)
	 */
	void LoadAliases(const char *aliasTypeList);
	void LoadPreferredAliases(const char *filename);

	void SummarizeModelCounts();

	/**
	 * Update the database with new variation filename (allow users to move the file to a permanent location specific to their system)
	 */
	void SetVariationFilename(const char *variationFilename);
	void ListGenes(std::ostream& os, Utility::StringArray& aliasList, Utility::StringArray& aliasType);
	/**
	 * Produce Gene/Gene models and write report to stream
	 * Returns the number of models produced
	 */
	void ProduceModels(std::ostream& os);

	void SetGeneExtension(uint geneBoundaryExt);
	void SetReportPrefix(const char *pref);
	void UseHtmlReports(bool useHtml);

	soci::session &DB();
	
	void LoadBuildConverter(const char* build);
	
	std::vector<uint> ManagerIDs();
	Knowledge::GroupManagerDB& GroupManager(uint idx);
protected:
	uint GetPopID(const char *pop);

	std::string dbFilename;										///< The name of the database file
	soci::session sociDB;										///< The database
	Knowledge::RegionManagerDB regions;						///< The genes
	std::map<uint, Knowledge::GroupManagerDB> groups;	///< The knowedge meta groups
	Knowledge::SnpDataset dataset;							///< the data associated with the user
	uint varVersion;												///< The variation version (to guarantee that the variations file is correct for the database being used)
	std::string variationFilename;							///< Filename for variation data-this might not be opened, depending on how the user loads their data
	std::stringstream reportLog;								///< The list of report filenames generated

	uint geneExtensionLength;									///< Length of extension on either side of a gene (not to be mixed with LD extension)
	bool htmlReports;												///< Turn on/off HTML report generation
	std::string reportPrefix;									///< Report Prefix
	std::multimap<uint, uint> geneLookup;
	Knowledge::GeneGeneModelArchive geneGeneModels;
public:
	static bool errorExit;										///< When exiting on errors, we won't report the files that "would" have been generated.
	LiftOver::ConverterDB buildConverter;					///< Converter structure for map files (using liftover chains)
};

inline
Application::Application() : dbFilename(""), varVersion(0), geneExtensionLength(0), htmlReports(false), reportPrefix("report") {}

inline
Knowledge::RegionManager *Application::GetRegions() {
	return &regions;
}

inline
Knowledge::SnpDataset *Application::GetDataset() {
	return &dataset;
}

inline
std::multimap<uint, uint> *Application::GetGeneLookup() {
	return &geneLookup;
}
inline
soci::session& Application::DB() {
	return sociDB;
}

inline
Knowledge::GeneGeneModelArchive *Application::GetGeneGeneModels() {
	return &geneGeneModels;
}
inline
void Application::SetGeneExtension(uint geneBoundaryExt) {
	geneExtensionLength = geneBoundaryExt;
}

inline
void Application::SetReportPrefix(const char *pref) {
	if (strcmp(pref, "") == 0)
		reportPrefix = "biofilter";
	else
		reportPrefix = pref;
}

inline
void Application::UseHtmlReports(bool doUse) {
	htmlReports = doUse;
}

inline
Application::~Application()  {
	if (!errorExit)
		std::cerr<<GetReportLog()<<"\n";
}

inline
std::vector<uint> Application::ManagerIDs() {
	std::vector<uint> ids;

	std::map<uint, Knowledge::GroupManagerDB>::iterator itr = groups.begin();
	std::map<uint, Knowledge::GroupManagerDB>::iterator end = groups.end();
	while (itr != end) 
		ids.push_back(itr++->first);
	return ids;
}

inline
Knowledge::GroupManagerDB& Application::GroupManager(uint idx) {
	return groups[idx];
}


}


#endif	/* APPLICATION_H */

