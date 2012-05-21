/*
 * InformationSQLite.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#include "InformationSQLite.h"

#include <sstream>
#include <stdlib.h>

using std::stringstream;

namespace Knowledge{

InformationSQLite::InformationSQLite(const string& filename) : _self_open(true){
	sqlite3_open(filename.c_str(), &_db);
}

InformationSQLite::InformationSQLite(sqlite3* db) : _db(db), _self_open(false){}

InformationSQLite::~InformationSQLite(){
	if (_self_open){
		sqlite3_close(_db);
	}
}

int InformationSQLite::getPopulationID(const string& pop_str){
	string queryStr = string("SELECT population_id FROM population "
			"WHERE population='") + pop_str + string("')");

	string result;
	int err_code = sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery,
			&result, NULL);

	if (err_code != 0){
		string queryStr = string("SELECT population_id FROM population "
				"WHERE population='n/a'");
		if(sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery, &result, NULL)){
			//NOTE: I should never get here in a properly formatted LOKI 2.0 database!
			return 1;
		}
	}
	return atoi(result.c_str());
}

const string InformationSQLite::getResourceVersion(const string& resource){
	string queryStr = string("SELECT version FROM versions WHERE element='") +
			resource + string("'");

	string result;
	int err_code = sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery,
			&result, NULL);

	if (err_code != 0){
		return "";
	}
	return result;
}

void InformationSQLite::getGroupTypes(const set<uint>& type_ids,
		map<int, string>& group_types_out){
	string query_str = string("SELECT group_type_id, group_type FROM group_type");

	if (type_ids.size() != 0){
		set<uint>::const_iterator itr = type_ids.begin();
		set<uint>::const_iterator end = type_ids.end();
		query_str += " WHERE group_type_id IN (";
		stringstream id_stream;
		id_stream << *itr;
		while(++itr != end){
			id_stream << "," << *itr;
		}
		query_str += id_stream.str() + ")";
	}

	sqlite3_exec(_db, query_str.c_str(), parseGroupTypeQuery, &group_types_out, NULL);
}

int InformationSQLite::parseSingleStringQuery(void* obj, int n_cols,
		char** col_vals, char** col_names){

	// Incorrect number of columns
	if (n_cols > 1){
		return 2;
	}

	// Not enough data!  Make sure to handle this case, because I think it
	// will actually come up.
	if (n_cols == 0){
		return 1;
	}
	string* ret_str_p = (string*) (obj);
	(*ret_str_p) = col_vals[0];
	return 0;
}

int InformationSQLite::parseGroupTypeQuery(void* obj, int n_cols,
		char** col_vals, char** col_names){

	if (n_cols != 2){
		return 2;
	}

	map<int, string>* map_obj = (map<int, string>*) obj;

	(*map_obj)[atoi(col_vals[0])] = col_vals[1];

	return 0;
}

}



