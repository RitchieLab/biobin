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
	GroupCollectionSQLite(RegionCollection& reg, const string& fn);
	GroupCollectionSQLite(RegionCollection& reg, sqlite3 *db_conn);

	virtual ~GroupCollectionSQLite();

	virtual void Load(const vector<string>& group_names,
			const unordered_set<uint>& ids);

protected:
	virtual uint getMaxGroup();

private:
	bool _self_open;
	// sqlite connection
	sqlite3 *_db;

	sqlite3_stmt* _group_name_stmt;

	Group* addGroup(sqlite3_stmt* group_query);

	void initQueries();

	static int parseGroupQuery(void*, int, char**, char**);
	static int parseGroupRelationshipQuery(void*, int, char**, char**);
	static int parseGroupAssociationQuery(void*, int, char**, char**);
	static int parseMaxGroupQuery(void*, int, char**, char**);
};

}




#endif /* GROUPCOLLECTIONSQLITE_H_ */
