/*
 * GroupCollectionSQLite.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include <sstream>
#include <stdlib.h>
#include <utility>

#include "GroupCollectionSQLite.h"

using std::stringstream;
using std::pair;

namespace Knowledge{

GroupCollectionSQLite::GroupCollectionSQLite(uint type, const string& name,
		const string& fn) :
			GroupCollection(type, name), _self_open(true){
	sqlite3_open(fn.c_str(),&_db);
}

GroupCollectionSQLite::GroupCollectionSQLite(uint type, const string& name,
		sqlite3 *db_conn) :
			GroupCollection(type, name), _self_open(false), _db(db_conn){}

GroupCollectionSQLite::~GroupCollectionSQLite(){
	if (_self_open){
		sqlite3_close(_db);
	}
}

uint GroupCollectionSQLite::Load(RegionCollection& regions, const unordered_set<uint>& ids){

	stringstream where_stream;
	where_stream << "WHERE groups.group_type_id = " << _id;

	if (ids.size() != 0){
		unordered_set<uint>::const_iterator itr = ids.begin();
		unordered_set<uint>::const_iterator end = ids.end();
		where_stream << " AND groups.group_id IN (" << *itr;
		while (++itr != end){
			where_stream << "," << *itr;
		}
		where_stream << ")";
	}

	string where_clause = where_stream.str();

	string group_query = "SELECT group_id, group_name, group_desc FROM groups " +
			where_clause;

	uint result = sqlite3_exec(_db, group_query.c_str(), parseGroupQuery, this, NULL);

	if (result != 0){
		string relationship_query = "SELECT parent_id, child_id FROM group_relationships "
				" INNER JOIN groups ON group_relationships.parent_id = groups.group_id " +
				where_clause;

		result = sqlite3_exec(_db, relationship_query.c_str(),
				parseGroupRelationshipQuery, this, NULL);
	}

	if (result != 0){
		string association_query = "SELECT group_id, gene_id FROM group_associations "
				"INNER JOIN groups USING (group_id) " + where_clause;

		pair<GroupCollection*, RegionCollection*> input =
				std::make_pair(this, &regions);

		result = sqlite3_exec(_db, association_query.c_str(),
				parseGroupAssociationQuery, &input, NULL);
	}

	return result;

}

int GroupCollectionSQLite::parseGroupQuery(void* obj, int n_cols, char** col_vals, char** col_names){

	GroupCollection *groups = (GroupCollection*) obj;

	// wrong # of columns!!
	if(n_cols != 3){
		return 2;
	}

	int group_id = atoi(col_vals[0]);
	groups->addGroup(group_id, col_vals[1], col_vals[2]);
	return 0;
}

int GroupCollectionSQLite::parseGroupRelationshipQuery(void* obj, int n_cols, char** col_vals, char** col_names){
	GroupCollection *groups = (GroupCollection*) obj;

	// wrong # of columns!!
	if(n_cols != 2){
		return 2;
	}

	int parent_id = atoi(col_vals[0]);
	int child_id = atoi(col_vals[1]);

	groups->addRelationship(parent_id, child_id);

	return 0;
}

int GroupCollectionSQLite::parseGroupAssociationQuery(void* obj, int n_cols, char** col_vals, char** col_names){
	pair<GroupCollection*, RegionCollection*>* input =
			(pair<GroupCollection*, RegionCollection*>*) obj;

	if(n_cols != 2){
		return 2;
	}

	int group_id = atoi(col_vals[0]);
	int region_id = atoi(col_vals[1]);

	input->first->addAssociation(group_id, region_id, *(input->second));

	return 0;
}

} //namespace Knowledge;
