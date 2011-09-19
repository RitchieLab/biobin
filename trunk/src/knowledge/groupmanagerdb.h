/* 
 * File:   groupmanagerdb.h
 * Author: torstees
 *
 * Created on March 23, 2011, 11:22 AM
 */

#ifndef GROUPMANAGERDB_H
#define	GROUPMANAGERDB_H

#include "groupmanager.h"
#include <soci.h>

namespace Knowledge {

class GroupManagerDB : public GroupManager {
public:
	GroupManagerDB() : GroupManager((uint)-1, MetaGroup::UNRECOGNIZED, "") {}
	GroupManagerDB(uint id, MetaGroup::Type groupType, const char *name = "??",
				const char *desc = "") : GroupManager(id, groupType, name, desc) {}
	GroupManagerDB(const GroupManagerDB& orig) : GroupManager(orig) {}
	virtual ~GroupManagerDB() {}

	/**
	 * Loads groups from the database whose IDs are found in the id list
    * @param sociDB
    * @param ids
    * @return
    */
	uint LoadFromDB(soci::session& sociDB,
			const Utility::IdCollection& ids,
			RegionManager& regions);
	/**
	 * Loads all regions from the database
	 */
	uint LoadFromDB(soci::session& sociDB, RegionManager& regions);

	/**
	 * Load the groups based on their name
    * @param sociDB Database connection
    * @param regions	The region information
    * @param groupNames The array of names to be loaded
	 * @param groupsFound Array of names that were found
    * @return Number of groups found
    */
	uint LoadFromDB(soci::session& sociDB,
	 const Utility::IdCollection& ids,
	 RegionManager& regions,
	 Utility::StringArray& groupNames);
private:

};




}

#endif	/* GROUPMANAGERDB_H */

