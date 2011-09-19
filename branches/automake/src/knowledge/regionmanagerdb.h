/* 
 * File:   dbregionmanager.h
 * Author: torstees
 *
 *
 * Created on March 9, 2011, 2:45 PM
 */

#ifndef DBREGIONMANAGER_H
#define	DBREGIONMANAGER_H

#include "regionmanager.h"
#include "snpdataset.h"
#include <soci.h>

namespace Knowledge {

class RegionManagerDB : public RegionManager {
public:
	RegionManagerDB() {}
	RegionManagerDB(const RegionManagerDB& orig) : RegionManager(*this) {}
	virtual ~RegionManagerDB() {}

	uint LoadFromDB(soci::session& sociDB, uint popID, Utility::StringArray& aliasList, Utility::IdCollection& ids, std::map<std::string, uint>& aliasToID);

	/**
	 * Loads regions from the database whose IDs are found in the id list
    * @param sociDB
    * @param ids
    * @return
    */
	uint LoadFromDB(soci::session& sociDB, uint popID, const Utility::IdCollection& ids );

	/**
	 * Loads all regions from the database
	 */
	uint LoadFromDB(soci::session& sociDB, uint popID);
	
	void GenerateLookupTable(std::map<uint, uint>& lookup);

	uint LoadFromDB(soci::session& sociDB, uint popID, Utility::StringArray& aliasList, std::map<std::string, uint>& aliasToID);
	

	/**
	 * Loads aliases for a given region. ID in this case is the internal geneID, not the index
	 */
	void LoadRegionAliases(soci::session& sociDB, Utility::StringArray& aliasList, std::map<std::string, uint>& aliasToId);

private:

};



}

#endif	/* DBREGIONMANAGER_H */

