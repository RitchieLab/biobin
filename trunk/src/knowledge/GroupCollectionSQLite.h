/*
 * GroupCollectionSQLite.h
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_GROUPCOLLECTIONSQLITE_H
#define KNOWLEDGE_GROUPCOLLECTIONSQLITE_H

#include "GroupCollection.h"

#include <sqlite3.h>

namespace Knowledge{

class GroupCollectionSQLite : public GroupCollection {

public:
	GroupCollectionSQLite(RegionCollection& reg, const std::string& fn, Information* info=0);
	GroupCollectionSQLite(RegionCollection& reg, sqlite3 *db_conn, Information* info=0);

	virtual ~GroupCollectionSQLite();

	virtual void Load(const std::vector<std::string>& group_names,
			const boost::unordered_set<unsigned int>& ids);

protected:
	virtual unsigned int getMaxGroup();

private:
	bool _self_open;
	bool _self_info;
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
