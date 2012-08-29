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
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

using std::vector;
using std::stringstream;
using std::ifstream;
using boost::is_any_of;
using boost::algorithm::split;

namespace Knowledge{

string InformationSQLite::_role_region_tbl = "__tmp_role_region";
string InformationSQLite::_role_zone_tbl = "__tmp_role_zone";

InformationSQLite::InformationSQLite(const string& filename) : _self_open(true){
	sqlite3_open(filename.c_str(), &_db);
	prepRoleTables();
	prepRoleStmt();
}

InformationSQLite::InformationSQLite(sqlite3* db) : _db(db), _self_open(false){
	prepRoleTables();
	prepRoleStmt();
}

InformationSQLite::~InformationSQLite(){
	sqlite3_finalize(_role_stmt);
	sqlite3_finalize(_region_role_stmt);

	if (_self_open){
		sqlite3_close(_db);
	}
}

int InformationSQLite::getPopulationID(const string& pop_str){
	string queryStr = string("SELECT ldprofile_id FROM ldprofile "
			"WHERE ldprofile='") + pop_str + string("'");

	string result;
	int err_code = sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery,
			&result, NULL);

	if (err_code != 0){
		string queryStr = string("SELECT ldprofile_id FROM ldprofile "
				"WHERE ldprofile=''");
		if(sqlite3_exec(_db, queryStr.c_str(), parseSingleStringQuery, &result, NULL)){
			//NOTE: I should never get here in a properly formatted LOKI 2.0 database!
			return 1;
		}
	}
	return atoi(result.c_str());
}

int InformationSQLite::getZoneSize(){

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

	// Look up dbSNP role here
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

	// Look up region role here
	int chr_idx = sqlite3_bind_parameter_index(_region_role_stmt, ":chr");
	int pos_idx = sqlite3_bind_parameter_index(_region_role_stmt, ":pos");
	int gid_idx = sqlite3_bind_parameter_index(_region_role_stmt, ":gid");
	sqlite3_bind_int(_region_role_stmt, chr_idx, loc.getChrom());
	sqlite3_bind_int(_region_role_stmt, pos_idx, loc.getPos());
	sqlite3_bind_int(_region_role_stmt, gid_idx, reg.getID());
	while(sqlite3_step(_region_role_stmt)==SQLITE_ROW){
		int role = sqlite3_column_int(_region_role_stmt, 0);
		db_role = _role_map.find(role);
		if (db_role != _role_map.end()){
			ret_val |= (*db_role).second;
		}
	}
	sqlite3_reset(_region_role_stmt);


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

void InformationSQLite::loadRoles(){
	// A mapping of roled to their ID
	map<string, int> role_id_map;
	map<string, int>::const_iterator role_itr;
	int max_role = 0;

	string max_role_sql = "SELECT max(role_id) FROM role";
	sqlite3_exec(_db, max_role_sql.c_str(), parseSingleIntQuery, &max_role, NULL);

	string insert_sql_null = "INSERT OR IGNORE INTO " + _role_region_tbl + " "
			"(role_id, chr, posMin, posMax) VALUES (?,?,?,?)";
	string insert_sql_gene = "INSERT OR IGNORE INTO " + _role_region_tbl + " "
			"(role_id, chr, posMin, posMax, biopolymer_id) VALUES (?,?,?,?,?)";

	string get_gene_sql = "SELECT biopolymer_id from biopolymer_name "
			"WHERE name=? AND source_id IN " + getSourceList();


	sqlite3_stmt* insert_stmt_gene;
	sqlite3_prepare_v2(_db, insert_sql_gene.c_str(), -1, &insert_stmt_gene, NULL);
	sqlite3_stmt* insert_stmt_null;
	sqlite3_prepare_v2(_db, insert_sql_null.c_str(), -1, &insert_stmt_null, NULL);

	sqlite3_stmt* find_gene_stmt;
	sqlite3_prepare_v2(_db, get_gene_sql.c_str(), -1, &find_gene_stmt, NULL);


	// Get rid of any indexes in the 2 tables of interest
	map<string, string> region_idx;
	getAndDropIndexes(_role_region_tbl, region_idx);

	// number of roles must be less than size of a long
	unsigned int n_roles = 3;

	vector<string>::const_iterator fn_itr = c_role_files.begin();
	while(fn_itr != c_role_files.end()){

		// for each role file, open it and load it into the db

		// Find the file and upload it...
		ifstream data_file((*fn_itr).c_str());
		if (!data_file.is_open()) {
			std::cerr << "WARNING: cannot find " << *fn_itr << ", ignoring.\n";
		} else {
			string line;
			vector<string> result;
			while (data_file.good()) {
				getline(data_file, line);
				split(result, line, is_any_of(" \n\t"), boost::token_compress_on);
				if((result.size() == 4 || result.size() == 5) && result[0][0] != '#'){
					short chr = Locus::getChrom(result[0]);
					int posMin = atoi(result[2].c_str());
					int posMax = atoi(result[3].c_str());

					int role_id;
					role_itr = role_id_map.find(result[1]);
					if(role_itr == role_id_map.end()){
						role_id = ++max_role;
						role_id_map[result[1]] = role_id;
						_role_map.insert(std::make_pair(role_id, getRole(result[1])));

					}else{
						role_id = (*role_itr).second;
					}

					if (result.size() == 4){
						if(chr != -1 && posMin && posMax){
							sqlite3_bind_int(insert_stmt_null, 1, role_id);
							sqlite3_bind_int(insert_stmt_null, 2, chr);
							sqlite3_bind_int(insert_stmt_null, 3, posMin);
							sqlite3_bind_int(insert_stmt_null, 4, posMax);

							while(sqlite3_step(insert_stmt_null)==SQLITE_ROW) {}

							sqlite3_reset(insert_stmt_null);
						}
					}else{
						// result size must be 5!
						sqlite3_bind_int(insert_stmt_gene, 1, role_id);
						sqlite3_bind_int(insert_stmt_gene, 2, chr);
						sqlite3_bind_int(insert_stmt_gene, 3, posMin);
						sqlite3_bind_int(insert_stmt_gene, 4, posMax);

						sqlite3_bind_text(find_gene_stmt, 1, result[4].c_str(), -1, SQLITE_STATIC);
						while(sqlite3_step(find_gene_stmt) == SQLITE_ROW){
							sqlite3_bind_int(insert_stmt_gene, 5, sqlite3_column_int(find_gene_stmt, 0));
							while(sqlite3_step(insert_stmt_gene) == SQLITE_ROW){}
							sqlite3_reset(insert_stmt_gene);
						}
						sqlite3_reset(find_gene_stmt);

					}
				}
			}
			data_file.close();
		}

		++fn_itr;

	}

	// If too many roles, print a warning
	if(n_roles > sizeof(unsigned long)){
		std::cerr << "WARNING: Too many roles.  "
				<< "Number of unique roles must be less than "
				<< sizeof(unsigned long) - 2 << "\n";
	}

	// Reload the indexes in the 2 tables of interest
	restoreIndexes(_role_region_tbl, region_idx);
	sqlite3_finalize(insert_stmt_gene);
	sqlite3_finalize(insert_stmt_null);
	sqlite3_finalize(find_gene_stmt);

	// Populate the zone table
	UpdateZones();

}

void InformationSQLite::UpdateZones(){
	// See loki_updater.py for the algorithm here
	// NOTE: blatantly and unabashedly copied from ldsplineimporter

	// Get the zone size
	int zone_size=getZoneSize();

	// Reverse any regions that are backwards
	string region_reverse_sql = "UPDATE " + _role_region_tbl +
			"SET posMin = posMax, posMax = posMin WHERE posMin > posMax";
	sqlite3_exec(_db, region_reverse_sql.c_str(), NULL, NULL, NULL);

	// Find the minimum and maximum zone sizes
	int minPos, maxPos = 0;
	string min_pos_sql = "SELECT MIN(posMin) FROM " + _role_region_tbl;
	string max_pos_sql = "SELECT MAX(posMax) FROM " + _role_region_tbl;

	sqlite3_exec(_db, min_pos_sql.c_str(), &parseSingleIntQuery, &minPos, NULL);
	sqlite3_exec(_db, max_pos_sql.c_str(), &parseSingleIntQuery, &maxPos, NULL);

	minPos = minPos / zone_size;
	maxPos = maxPos / zone_size;

	// Insert all zones needed into a temporary table
	string tmp_zone_tbl = "__zone_tmp";
	string zone_tbl_sql = "CREATE TEMPORARY TABLE " + tmp_zone_tbl + " "
			"(zone INTEGER PRIMARY KEY NOT NULL)";
	sqlite3_exec(_db, zone_tbl_sql.c_str(), NULL, NULL, NULL);

	string zone_insert_query = "INSERT INTO " + tmp_zone_tbl + " VALUES (?)";
	sqlite3_stmt* zone_ins_stmt;
	sqlite3_prepare_v2(_db, zone_insert_query.c_str(), -1, &zone_ins_stmt, NULL);

	for (int i=minPos; i <= maxPos; i++){
		sqlite3_bind_int(zone_ins_stmt, 1, i);
		while(sqlite3_step(zone_ins_stmt) == SQLITE_ROW){}
		sqlite3_reset(zone_ins_stmt);
	}
	sqlite3_finalize(zone_ins_stmt);

	// drop the indexes on zone table
	map<string, string> index_map;
	getAndDropIndexes(_role_zone_tbl,index_map);

	// Get rid of the entire region zone table
	string zone_del_sql = "DELETE FROM " + _role_zone_tbl;
	sqlite3_exec(_db, zone_del_sql.c_str(), NULL, NULL, NULL);

	// OK, now we're ready to do some insertin'!
	stringstream zone_ins_str;
	zone_ins_str << "INSERT OR IGNORE INTO " + _role_zone_tbl + " (role_region_id,chr,zone) "
			<< "SELECT rb.role_region_id, rb.chr, z.zone "
			<< "FROM " << _role_region_tbl << " AS rb JOIN " << tmp_zone_tbl << " AS z "
			<< "ON z.zone >= rb.posMin / " << zone_size << " "
			<< "AND z.zone <= rb.posMax / " << zone_size << " ";

	string zone_ins_sql = zone_ins_str.str();
	sqlite3_exec(_db, zone_ins_sql.c_str(), NULL, NULL, NULL);

	string tmp_tbl_drop = "DROP TABLE " + tmp_zone_tbl;
	sqlite3_exec(_db, tmp_tbl_drop.c_str(), NULL, NULL, NULL);

	restoreIndexes(_role_zone_tbl, index_map);
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
			"('utr-3','utr-5','regulatory')";
	sqlite3_exec(_db, reg_sql.c_str(), parseMultiIntQuery, &role_ids, NULL);

	role_itr = role_ids.begin();
	while(role_itr != role_ids.end()){
		_role_map.insert(std::make_pair(*role_itr,REGULATORY));
		++role_itr;
	}
	role_ids.clear();

	// Prep the SQL statement to get the code
	string role_sql = "SELECT role_id FROM snp_locus "
			"INNER JOIN snp_biopolymer_role USING (rs) "
			"WHERE chr=? AND pos=? AND biopolymer_id=? AND "
			"snp_biopolymer_role.source_id IN " + getSourceList();
	sqlite3_prepare_v2(_db, role_sql.c_str(), -1, &_role_stmt, NULL);

	// Prep the SQL statement for the region
	stringstream region_role_sql_str;
	region_role_sql_str << "SELECT role_id FROM " << _role_zone_tbl << " AS z "
			<< "INNER JOIN " << _role_region_tbl << " AS r USING (role_region_id) "
			<< "WHERE z.chr=:chr AND z.zone=:pos/" << getZoneSize() << " "
			<< "AND posMin<=:pos AND posMax>=:pos AND "
			<< "(biopolymer_id IS NULL OR biopolymer_id=:gid)";

	sqlite3_prepare_v2(_db, region_role_sql_str.str().c_str(), -1, &_region_role_stmt, NULL);

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
	if(col_vals[0]){
		int* ret_int_p = (int*) (obj);
		(*ret_int_p) = atoi(col_vals[0]);
	}
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

void InformationSQLite::getAndDropIndexes(const string& tbl_name,
		map<string, string>& indexes_out) {

	string idx_cmd = "SELECT name, sql FROM sqlite_master "
			"WHERE type='index' AND tbl_name='" + tbl_name + "' "
					"AND sql NOT NULL";

	sqlite3_exec(_db, idx_cmd.c_str(), &parseRegionIndex, &indexes_out, NULL);

	// Now, drop those indexes!const string& tbl_name,
	string drop_cmd = "DROP INDEX ";

	map<string, string>::const_iterator idx_itr = indexes_out.begin();
	while (idx_itr != indexes_out.end()) {
		string drop_tbl = "'" + (*idx_itr).first + "'";
		string sql_str = drop_cmd + drop_tbl;
		sqlite3_exec(_db, (drop_cmd + drop_tbl).c_str(), NULL, NULL, NULL);
		++idx_itr;
	}
}

void InformationSQLite::restoreIndexes(const string& tbl_name, const map<string, string>& index_map){
	// Recretate the indexes
	map<string, string>::const_iterator idx_itr = index_map.begin();
	while(idx_itr != index_map.end()){
		sqlite3_exec(_db, (*idx_itr).second.c_str(), NULL, NULL, NULL);
		++idx_itr;
	}
	//Analyze the table
	string analyze_sql = "ANALYZE '" + tbl_name + "'";
	sqlite3_exec(_db, analyze_sql.c_str(), NULL, NULL, NULL);

}

void InformationSQLite::prepRoleTables(){

	vector<string> sql_stmts;
	sql_stmts.push_back("CREATE TEMPORARY TABLE " + _role_region_tbl +
			"("
			"role_region_id INTEGER PRIMARY KEY AUTOINCREMENT,"
			"role_id TINYINT NOT NULL,"
			"chr TINYINT NOT NULL,"
			"posMin BIGINT NOT NULL,"
			"posMax BIGINT NOT NULL,"
			"biopolymer_id INTEGER"
			")");
	sql_stmts.push_back("CREATE TEMPORARY TABLE " + _role_zone_tbl +
			"("
			"role_region_id INTEGER NOT NULL,"
			"chr TINYINT NOT NULL,"
			"zone INTEGER NOT NULL,"
			"PRIMARY KEY (role_region_id,chr,zone)"
			")");

	sql_stmts.push_back("CREATE INDEX '" + _role_region_tbl + "__chr_min' ON " + _role_region_tbl + " (chr,posMin)");
	sql_stmts.push_back("CREATE INDEX '" + _role_region_tbl + "__chr_max' ON " + _role_region_tbl + " (chr,posMax)");
	sql_stmts.push_back("CREATE INDEX '" + _role_zone_tbl + "__zone' ON " + _role_zone_tbl + " (chr,zone,role_region_id)");

	vector<string>::const_iterator sql_itr = sql_stmts.begin();
	while(sql_itr != sql_stmts.end()){
		sqlite3_exec(_db, (*sql_itr).c_str(), NULL, NULL, NULL);
		++sql_itr;
	}

}

int InformationSQLite::parseRegionIndex(void* obj, int n_cols, char** col_vals, char** col_names){
	if (n_cols != 2){
		return 2;
	}

	map<string, string>* idx_map_p = static_cast<map<string, string>*>(obj);
	(*idx_map_p)[col_vals[0]] = col_vals[1];
	return 0;
}

}



