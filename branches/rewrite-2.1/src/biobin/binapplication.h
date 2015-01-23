/* 
 * File:   binapplication.h
 * Author: torstees
 *
 * Created on June 22, 2011, 10:35 AM
 * 
 * One of the big TODO list things would be to integrate the two 
 * forms of SNPs: Biofilter and BioBin. Right now, we have two
 * very different approaches to SNPs. For the biofilter, we only
 * need a way to recognize names and associate them with a basepair
 * and chromosome. However, for biobin, we need to maintain alleles
 * and provide the ability to perform genotype conversion and 
 * some other stuff. So, the Locus object is much more complex. 
 * 
 * Most likely it's just a matter of moving the biobin Locus class
 * to someplace common and changing the biofilter code to use it...
 * 
 */

#ifndef BINAPPLICATION_H
#define	BINAPPLICATION_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <new>
#include <deque>
#include <utility>
#include <set>
#include <cstdio>

#include <sqlite3.h>

// Use the boost filesystem library to work with OS-independent paths
#include <boost/filesystem.hpp>

#include "Bin.h"
#include "binmanager.h"
#include "PopulationManager.h"

#include "knowledge/Locus.h"
#include "knowledge/Region.h"
#include "knowledge/Group.h"
#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Information.h"
#include "knowledge/liftover/Converter.h"
#include "knowledge/GroupCollectionSQLite.h"
#include "knowledge/liftover/ConverterSQLite.h"

namespace BioBin {
	
class BinApplication{
public:
	BinApplication(const std::string& db_fn, const std::string& vcf_file);
	~BinApplication();

	// We want to separate the next two steps because some jobs require only one
	// of the two be performed, and they do take time to complete
	/**
	 * Loads region data from database
    * @param pop indicates LD boundaries to be used, if any (if NO-LD, then the variable geneExtensionLength is used for static extensions
    * @param regionAliases A list of region aliases to be loaded. This will restrict the genes loaded to only those found (this is intended for certain types of reporting only)
    * @return Number of regions loaded
    */
	template <class T_cont>
	unsigned int LoadRegionData(T_cont& aliasesNotFound,
			const std::vector<std::string>& aliasList);		///< Comma separated list of strings

	/**
	 * Loads group data using group name rather than IDs
    * @param userDefinedGroups
    * @param groupNames
    * @return
    */
	template <class T1_cont>
	unsigned int LoadGroupDataByName(T1_cont& userDefinedGroups);

	unsigned int GetPopulationID(const std::string& pop){return _info->getPopulationID(pop);}

	void SetGeneExtension(unsigned int geneBoundaryExt){geneExtensionLength = geneBoundaryExt;}
	void SetReportPrefix(const std::string& pref){reportPrefix=(pref=="")?"biobin":pref;}
	void InitVcfDataset(const std::string& genomicBuild);
	void loadRoles(){_info->loadRoles(*regions);}

	/**
	 * Initialize the bins.  After this call, the binmanager will have the final
	 * bins set up and ready to output.
     */
	void InitBins();

	void writeBinData(const std::string& filename, const std::string& sep=",") const;
	void writeLoci(const std::string& filename, const std::string& sep=",") const;

	static bool c_transpose_bins;
	static bool errorExit;										///< When exiting on errors, we won't report the files that "would" have been generated.
	static std::string reportPrefix;
	static bool c_print_sources;
	static bool c_print_populations;
	static bool s_run_normal;

private:
	void Init(const std::string& dbFilename, bool reportVersions);

	void printEscapedString(std::ostream& os, const std::string& toPrint, const std::string& toRepl, const std::string& replStr) const;
	std::string getEscapeString(const std::string& sep) const;

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
	std::deque<Knowledge::Locus*> dataset;

	std::vector<FILE*> _locus_bins;

	///< The variation version (to guarantee that the variations file is correct for the database being used)
	unsigned int varVersion;
	///< Filename for variation data-this might not be opened, depending on how the user loads their data
	std::string variationFilename;

	///< Length of extension on either side of a gene (not to be mixed with LD extension)
	unsigned int geneExtensionLength;

	PopulationManager _pop_mgr;

	BinManager* binData;

//Everything from here on down has to do with installing a new handler that
// will try to get sqlite to give up some of its cache
public:
	// A function to release as much db cache as possible
	static void releaseDBCache();

private:

	static std::new_handler currentHandler;

private:


};


template <class T_cont>
unsigned int BinApplication::LoadRegionData(T_cont& aliasesNotFound, const std::vector<std::string>& aliasList) {

	regions->Load(aliasList);

	//T_cont::const_iterator pos = aliasesNotFound.end();

	std::vector<std::string>::const_iterator itr = aliasList.begin();
	std::vector<std::string>::const_iterator end = aliasList.end();
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
unsigned int BinApplication::LoadGroupDataByName(T1_cont& userDefinedGroups) {

	groups = new Knowledge::GroupCollectionSQLite(*regions, _db, _info);
	groups->Load();

	std::vector<std::string> unmatchedAliases;
	typename T1_cont::const_iterator udItr = userDefinedGroups.begin();
	typename T1_cont::const_iterator udEnd = userDefinedGroups.end();

	while (udItr != udEnd) {
		groups->LoadArchive(*udItr, unmatchedAliases);
		++udItr;
	}

	return groups->size();

}

}
#endif	/* BINAPPLICATION_H */

