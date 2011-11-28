/*
 * RegionCollectionSQLite.h
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_REGIONCOLLECTIONSQLITE_H
#define KNOWLEDGE_REGIONCOLLECTIONSQLITE_H_

// Use the straight-up sqlite interface
#include <sqlite3.h>

#include "RegionCollection.h"


namespace Knowledge{

/**
 * A class that implements the region collection from a SQLite database
 */
class RegionCollectionSQLite : public RegionCollection{

public:
	RegionCollectionSQLite(const string& fn);
	RegionCollectionSQLite(sqlite3 *);
	virtual ~RegionCollectionSQLite();

	/**
	 * Loading function - must be subclassed
	 */
	virtual uint Load(const uint popID,
			const unordered_set<uint>& ids,
			const vector<string>& aliasList);

private:
	bool self_open;
	// sqlite connection
	sqlite3 *db;

	// callback functions for sqlite interface
	// NOTE: the 1st argument will be a pointer to a RegionCollection object
	static int parseRegionQuery(void*, int, char**, char**);
};

}




#endif /* REGIONCOLLECTIONSQLITE_H_ */
