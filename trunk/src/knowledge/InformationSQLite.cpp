/*
 * InformationSQLite.cpp
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#include "InformationSQLite.h"
#include "Locus.h"
#include "Region.h"

#include <sstream>
#include <stdlib.h>
#include <vector>

using std::vector;
using std::stringstream;

namespace Knowledge{

InformationSQLite::InformationSQLite(const string& filename) : _self_open(true){
	sqlite3_open(filename.c_str(), &_db);
	prepRoleStmt();
}

InformationSQLite::InformationSQLite(sqlite3* db) : _db(db), _self_open(false){
	prepRoleStmt();
}

InformationSQLite::~InformationSQLite(){
	sqlite3_finalize(_role_stmt);

	if (_self_open){
		sqlite3_close(_db);
	}
}

int InformationSQLite::getPopulationID(const string& pop_str){
	string queryStr = string("SELECT ldprofile_id FROM ldprofile "
			"WHERE ldprofile='") + pop_str + string("')");

	string result;
	int err_code = sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery,
			&result, NULL);

	if (err_code != 0){
		string queryStr = string("SELECT ldprofile_id FROM ldprofile "
				"WHERE ldprofile='n/a'");
		if(sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery, &result, NULL)){
			//NOTE: I should never get here in a properly formatted LOKI 2.0 database!
			return 1;
		}
	}
	return atoi(result.c_str());
}

int InformationSQLite::getZoneSize(){
	// TODO: Add db query to get the zone size

	// default is 100K
	int zone_size = 100000;
	string zone_sql = "SELECT value FROM setting "
			"WHERE setting='zone_size'";
	sqlite3_exec(_db, zone_sql.c_str(), parseSingleIntQuery, &zone_size, NULL);

	return zone_size;
}

void InformationSQLite::getGroupTypes(const set<uint>& type_ids,
		map<int, string>& group_types_out){
	string query_str = string("SELECT source_id, source FROM source");

	if (type_ids.size() != 0){
		set<uint>::const_iterator itr = type_ids.begin();
		set<uint>::const_iterator end = type_ids.end();
		query_str += " WHERE source_id IN (";
		stringstream id_stream;
		id_stream << *itr;
		while(++itr != end){
			id_stream << "," << *itr;
		}
		query_str += id_stream.str() + ")";
	}

	sqlite3_exec(_db, query_str.c_str(), parseGroupTypeQuery, &group_types_out, NULL);
}

int InformationSQLite::getSNPRole(const Locus& loc, const Region& reg){
	int ret_val = 0;
	map<int, Information::snp_role>::const_iterator db_role = _role_map.end();

	// look it up here
	sqlite3_bind_int(_role_stmt, 1, loc.getChrom());
	sqlite3_bind_int(_role_stmt, 2, loc.getPos());
	sqlite3_bind_int(_role_stmt, 3, reg.getID());
	while(sqlite3_step(_role_stmt)==SQLITE_ROW){
		int role = sqlite3_column_int(_role_stmt, 0);
		db_role = _role_map.find(role);
		if (db_role != _role_map.end()){
			ret_val |= (*db_role).second;
		}
	}
	sqlite3_reset(_role_stmt);

	return ret_val;
}

void InformationSQLite::printPopulations(ostream& os){
	string pop_sql = "SELECT ldprofile, comment FROM ldprofile";

	os << "Population\tComment\n";
	sqlite3_exec(_db, pop_sql.c_str(), printQueryResult, &os, NULL);
}

void InformationSQLite::printSources(ostream& os){
	string src_sql = "SELECT source_id, source FROM source";

	os << "ID\tSource Name\n";
	sqlite3_exec(_db, src_sql.c_str(), printQueryResult, &os, NULL);

}

void InformationSQLite::prepRoleStmt(){
	vector<int> role_ids;
	vector<int>::const_iterator role_itr;
	// intron codes
	string intron_sql = "SELECT role_id FROM role WHERE role IN "
			"('intron','splice-3','splice-5')";
	sqlite3_exec(_db, intron_sql.c_str(), parseMultiIntQuery, &role_ids, NULL);

	role_itr = role_ids.begin();
	while(role_itr != role_ids.end()){
		_role_map.insert(std::make_pair(*role_itr, INTRON));
		++role_itr;
	}
	role_ids.clear();

	// exon codes
	string exon_sql = "SELECT role_id FROM role WHERE role IN "
			"('cds-synon','stop-gain','missense','frameshift')";
	sqlite3_exec(_db, exon_sql.c_str(), parseMultiIntQuery, &role_ids, NULL);

	role_itr = role_ids.begin();
	while(role_itr != role_ids.end()){
		_role_map.insert(std::make_pair(*role_itr, EXON));
		++role_itr;
	}
	role_ids.clear();

	// regulatory codes
	string reg_sql = "SELECT role_id FROM role WHERE role IN "
			"('utr-3','utr-5')";
	sqlite3_exec(_db, reg_sql.c_str(), parseMultiIntQuery, &role_ids, NULL);

	role_itr = role_ids.begin();
	while(role_itr != role_ids.end()){
		_role_map.insert(std::make_pair(*role_itr,REGULATORY));
		++role_itr;
	}
	role_ids.clear();

	// Prep the SQL statement to get the code
	string role_sql = "SELECT role_id FROM snp_locus INNER JOIN snp_biopolymer_role USING (rs) "
			"WHERE chr=? AND pos=? AND biopolymer_id=?";
	sqlite3_prepare_v2(_db, role_sql.c_str(), -1, &_role_stmt, NULL);

}

const set<unsigned int>& InformationSQLite::getSourceIds(){
	if(_s_source_ids.size() == 0){

		vector<int> query_results;
		string source_query = "SELECT source_id FROM source ";
		string where_clause = "";

		if(c_source_names.size() != 0){
			stringstream where_str;
			where_str << "WHERE source IN (";
			for(unsigned int i=0; i<c_source_names.size(); i++){
				if(i){
					where_str << ",";
				}
				where_str << "'" << c_source_names[i] << "'";
			}
			where_str << ")";

			where_clause = where_str.str();
		}

		string source_sql = source_query + where_clause;
		sqlite3_exec(_db, source_sql.c_str(), parseMultiIntQuery, &query_results, NULL);

		_s_source_ids.insert(query_results.begin(), query_results.end());
	}
	return _s_source_ids;
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

int InformationSQLite::parseSingleIntQuery(void* obj, int n_cols,
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
	int* ret_int_p = (int*) (obj);
	(*ret_int_p) = atoi(col_vals[0]);
	return 0;
}

int InformationSQLite::parseMultiIntQuery(void* obj, int n_cols, char** col_vals, char** col_names){
	if (n_cols != 1){
		return 2;
	}
	vector<int>* int_obj = static_cast<vector<int>*>(obj);

	int_obj->push_back(atoi(col_vals[0]));

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

int InformationSQLite::printQueryResult(void* obj, int n_cols, char** col_vals, char** col_names){

	ostream & os = *static_cast<ostream *>(obj);
	os << col_vals[0];
	for (int i = 1; i < n_cols; i++){
		os << "\t" << col_vals[i];
	}
	os << "\n";

	return 0;
}

}



