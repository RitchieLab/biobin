/*
 * GroupCollectionSQLite.h
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_GROUPCOLLECTIONSQLITE_H
#define KNOWLEDGE_GROUPCOLLECTIONSQLITE_H

#include <sqlite3.h>

#include <string>

#include <boost/unordered_set.hpp>

#include "GroupCollection.h"

using std::string;
using boost::unordered_set;

namespace Knowledge{

class GroupCollectionSQLite : public GroupCollection {

public:
	GroupCollectionSQLite(uint type, const string& name, const string& fn);
	GroupCollectionSQLite(uint type, const string& name, sqlite3 *db_conn);

	virtual ~GroupCollectionSQLite();

	virtual uint Load(RegionCollection& regions, const unordered_set<uint>& ids);

private:
	bool _self_open;
	// sqlite connection
	sqlite3 *_db;

	static int parseGroupQuery(void*, int, char**, char**);
	static int parseGroupRelationshipQuery(void*, int, char**, char**);
	static int parseGroupAssociationQuery(void*, int, char**, char**);
};

}




#endif /* GROUPCOLLECTIONSQLITE_H_ */
