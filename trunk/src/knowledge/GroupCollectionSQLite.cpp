/*
 * GroupCollectionSQLite.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include <sstream>
#include <stdlib.h>
#include <utility>
#include <deque>

#include "GroupCollectionSQLite.h"
#include "RegionCollection.h"
#include "InformationSQLite.h"

using std::stringstream;
using std::pair;
using std::deque;

namespace Knowledge{

GroupCollectionSQLite::GroupCollectionSQLite(RegionCollection& reg,
		const string& fn, Information* info) :
			GroupCollection(reg), _self_open(true), _self_info(!info){
	sqlite3_open(fn.c_str(),&_db);
	getMaxGroup();
	if(!info){
		_info = new InformationSQLite(_db);
	}else{
		_info = info;
	}

	initQueries();
}

GroupCollectionSQLite::GroupCollectionSQLite(RegionCollection& reg,
		sqlite3 *db_conn, Information* info) :
			GroupCollection(reg), _self_open(false), _self_info(!info), _db(db_conn) {
	getMaxGroup();
	if(!info){
		_info = new InformationSQLite(_db);
	}else{
		_info = info;
	}

	initQueries();
}

GroupCollectionSQLite::~GroupCollectionSQLite(){
	sqlite3_finalize(_group_name_stmt);

	if (_self_open){
		sqlite3_close(_db);
	}
	if(_self_info){
		delete _info;
	}
}

void GroupCollectionSQLite::Load(const vector<string>& group_names,
		const unordered_set<uint>& ids){

	//First, get a list of IDs corresponding to the group names
	unordered_set<uint> id_list(ids);

	stringstream src_str;
	const set<unsigned int>& src_list = _info->getSourceIds();
	if (src_list.size() > 0){
		src_str << " AND 'group'.source_id IN (";
		set<unsigned int>::const_iterator src_itr = src_list.begin();
		src_str << *src_itr;
		while(++src_itr != src_list.end()){
			src_str << "," << *src_itr;
		}
		src_str << ") ";
	}

	if(group_names.size() > 0){


		string command = "SELECT group_id FROM group_name "
				"INNER JOIN group USING (group_id) "
				"WHERE name=:label" + src_str.str();
		sqlite3_stmt* name_stmt;
		sqlite3_prepare_v2(_db, command.c_str(), -1, &name_stmt, NULL);
		int label_idx = sqlite3_bind_parameter_index(name_stmt, ":label");

		vector<string>::const_iterator n_itr = group_names.begin();
		while(n_itr != group_names.end()){
			sqlite3_bind_text(name_stmt, label_idx, (*n_itr).c_str(), -1, SQLITE_STATIC);
			while(sqlite3_step(name_stmt)==SQLITE_ROW){
				id_list.insert(static_cast<uint>(sqlite3_column_int(name_stmt, 0)));
			}
			sqlite3_reset(name_stmt);
			++n_itr;
		}

		sqlite3_finalize(name_stmt);
	}


	vector<string>::const_iterator rel_type = child_types.begin();
	stringstream rel_stream;
	rel_stream << "('" << *rel_type << "'";
	while(++rel_type != child_types.end()){
		rel_stream << ", '" << *rel_type << "'";
	}
	rel_stream << ")";

	// OK, now get a list of IDs that are children of anything in the group,
	// but only if we are filtering
	if (id_list.size() > 0){

		// We're going to check all of the groups we have, and add the
		// new groups to the list here
		deque<uint> loadables(id_list.begin(), id_list.end());

		string child_cmd = "SELECT related_group_id FROM group_group "
				"INNER JOIN relationship USING (relationship_id) "
				"WHERE relationship IN " + rel_stream.str() + " "
				"AND direction=1 AND group_id=:gid";

		sqlite3_stmt* child_stmt;
		sqlite3_prepare_v2(_db, child_cmd.c_str(), -1, &child_stmt, NULL);
		int group_idx = sqlite3_bind_parameter_index(child_stmt, ":gid");

		while(loadables.size() > 0){
			sqlite3_bind_int(child_stmt, group_idx, loadables.front());
			while(sqlite3_step(child_stmt) == SQLITE_ROW){
				uint new_id = static_cast<uint>(sqlite3_column_int(child_stmt, 0));
				if (id_list.find(new_id) != id_list.end()){
					id_list.insert(new_id);
					loadables.push_back(new_id);
				}
			}
			sqlite3_reset(child_stmt);
			loadables.pop_front();
		}

		sqlite3_finalize(child_stmt);
	}

	string group_cmd = "SELECT group_id, label, description, source "
			"FROM group_biopolymer "
			"INNER JOIN 'group' USING (group_id) "
			"INNER JOIN source ON 'group'.source_id=source.source_id "
			"WHERE group_biopolymer.biopolymer_id=:region_id " + src_str.str();

	sqlite3_stmt* group_stmt;
	sqlite3_prepare_v2(_db, group_cmd.c_str(), -1, &group_stmt, NULL);
	int region_idx = sqlite3_bind_parameter_index(group_stmt, ":region_id");

	deque<uint> child_groups;

	RegionCollection::const_iterator r_itr = _regions.begin();
	while(r_itr != _regions.end()){
		sqlite3_bind_int(group_stmt, region_idx, (*r_itr)->getID());
		while(sqlite3_step(group_stmt)==SQLITE_ROW){
			uint group_id = static_cast<uint>(sqlite3_column_int(group_stmt, 0));
			if(id_list.size() == 0 || id_list.find(group_id) != id_list.end()){
				Group* gp;
				unordered_map<uint, Group*>::const_iterator test_itr = _group_map.find(group_id);
				unordered_map<uint, Group*>::const_iterator test_end = _group_map.end();
				if (_group_map.find(group_id) == _group_map.end()){
					child_groups.push_back(group_id);
					gp = addGroup(group_stmt);
				} else {
					gp = _group_map[group_id];
				}
				gp->addRegion(**r_itr);
				(*r_itr)->addGroup(*gp);
				_group_associations[group_id].insert(*r_itr);
			}
		}
		sqlite3_reset(group_stmt);
		++r_itr;
	}

	sqlite3_finalize(group_stmt);

	// OK, now load parents of all of the groups that contain the regions
	string parent_cmd = "SELECT group_group.group_id, 'group'.label, "
			"'group'.description, source "
			"FROM group_group INNER JOIN 'group' USING (group_id) "
			"INNER JOIN source ON 'group'.source_id='group'.source_id "
			"INNER JOIN relationship ON group_group.relationship_id=relationship.relationship_id "
			"WHERE group_group.related_group_id=:gid AND direction=-1 "
			"AND relationship IN " + rel_stream.str();

	sqlite3_stmt* parent_stmt;
	sqlite3_prepare_v2(_db, parent_cmd.c_str(), -1, &parent_stmt, NULL);
	int gp_idx = sqlite3_bind_parameter_index(parent_stmt, ":gid");

	while(child_groups.size() > 0){
		sqlite3_bind_int(parent_stmt, gp_idx, child_groups.front());
		while(sqlite3_step(parent_stmt)==SQLITE_ROW){
			uint parent_id = static_cast<uint>(sqlite3_column_int(parent_stmt,0));
			if(id_list.size() == 0 || id_list.find(parent_id) != id_list.end()){
				Group* parent;
				if (_group_map.find(parent_id) == _group_map.end()){
					parent = addGroup(parent_stmt);
					child_groups.push_back(parent_id);
				}else{
					parent = _group_map[parent_id];
				}
				parent->addChild(*_group_map[child_groups.front()]);
				_group_relationships[parent_id].insert(child_groups.front());
			}
		}
		sqlite3_reset(parent_stmt);
		child_groups.pop_front();
	}

	sqlite3_finalize(parent_stmt);
}

uint GroupCollectionSQLite::getMaxGroup() {
	if (_max_group == 0){
		string query_str = "SELECT IFNULL(MAX(rowid), 1) FROM 'group'";
		sqlite3_exec(_db, query_str.c_str(), parseMaxGroupQuery, &_max_group, NULL);
	}

	return _max_group;

}

Group* GroupCollectionSQLite::addGroup(sqlite3_stmt* group_query){

	uint group_id = static_cast<uint>(sqlite3_column_int(group_query, 0));
	
	// Get all of the aliases
	sqlite3_bind_int(_group_name_stmt, 1, group_id);
	stringstream alias_str;
	int n_alias=-1;
	while(sqlite3_step(_group_name_stmt) == SQLITE_ROW){
		if(++n_alias){
			alias_str << ",";
		}
		alias_str << sqlite3_column_text(group_query, 0);
	}
	sqlite3_reset(_group_name_stmt);

	string src_name = (const char *) (sqlite3_column_text(group_query, 3));
	string gp_name = src_name + ":" + (const char *) (sqlite3_column_text(group_query, 1));

	string gp_desc = "";
	if(sqlite3_column_type(group_query,2) != SQLITE_NULL){
		gp_desc = (const char *) (sqlite3_column_text(group_query, 2));
	}
	string gp_alias = (const char *) (sqlite3_column_text(group_query, 3));

	return AddGroup(group_id, gp_name, gp_desc);
}

void GroupCollectionSQLite::initQueries(){
	string grp_alias_cmd = "SELECT name from group_name WHERE group_id=?";
	sqlite3_prepare_v2(_db, grp_alias_cmd.c_str(), -1, &_group_name_stmt, NULL);
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
	input->first->addAssociation(group_id, region_id);

	return 0;
}

int GroupCollectionSQLite::parseMaxGroupQuery(void* obj, int n_cols,
		char** col_vals, char** col_names){

	if (n_cols != 1){
		return 2;
	}

	int* result = (int*) obj;
	(*result) = atoi(col_vals[0]);
	return 0;
}

} //namespace Knowledge;
