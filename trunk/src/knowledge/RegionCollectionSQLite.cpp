/*
 * RegionCollectionSQLite.cpp
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#include "RegionCollectionSQLite.h"
#include "Locus.h"

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>


#include "Locus.h"

using std::stringstream;
using std::ifstream;
using boost::algorithm::split;
using boost::algorithm::is_any_of;

namespace Knowledge{

string RegionCollectionSQLite::_s_tmp_region_tbl = "__tmp_region";
string RegionCollectionSQLite::_s_tmp_zone_tbl = "__tmp_zone";


RegionCollectionSQLite::~RegionCollectionSQLite(){
	sqlite3_finalize(_region_name_stmt);
	sqlite3_finalize(_region_bound_stmt);

	if (self_open){
		sqlite3_close(db);
	}
}

void RegionCollectionSQLite::loadFiles(){
	// Only do this if we have something to load
	if (c_region_files.size() == 0){
		return;
	}

	string tmp_region_sql = "CREATE TEMPORARY TABLE " + _s_tmp_region_tbl + " ("
			"region_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
			"label VARCHAR(64) NOT NULL,"
			"chr TINYINT NOT NULL,"
			"posMin BIGINT NOT NULL,"
			"posMax BIGINT NOT NULL"
		")";

	sqlite3_exec(db, tmp_region_sql.c_str(), NULL, NULL, NULL);

	string tmp_zone_sql = "CREATE TEMPORARY TABLE " + _s_tmp_zone_tbl + " ("
			"region_id INTEGER NOT NULL,"
			"chr TINYINT NOT NULL,"
			"zone INTEGER NOT NULL,"
			"PRIMARY_KEY (region_id, chr, zone)"
		")";

	sqlite3_exec(db, tmp_zone_sql.c_str(), NULL, NULL, NULL);

	// Insert a phony record to set the autoincrement
	string region_id_sql = "INSERT INTO " + _s_tmp_region_tbl + " "
			"(region_id, label, chr, posMin, posMax) "
			"SELECT MAX(biopolymer_id), 'Fake', -1, 0, 0 FROM biopolymer";

	sqlite3_exec(db, region_id_sql.c_str(), NULL, NULL, NULL);

	vector<string>::const_iterator fn_itr = c_region_files.begin();
	while(fn_itr != c_region_files.end()){
		loadFile(*fn_itr);
		++fn_itr;
	}

	updateZones();

	string reg_min_index_sql = "CREATE INDEX " + _s_tmp_region_tbl + "_posMin "
			"ON " + _s_tmp_region_tbl + " (chr, posMin)";
	string reg_max_index_sql = "CREATE INDEX " + _s_tmp_region_tbl + "_posMax "
				"ON " + _s_tmp_region_tbl + " (chr, posMax)";

	sqlite3_exec(db, reg_min_index_sql.c_str(), NULL, NULL, NULL);
	sqlite3_exec(db, reg_max_index_sql.c_str(), NULL, NULL, NULL);

}

void RegionCollectionSQLite::loadFile(const string& fn){

	string insert_sql = "INSERT OR IGNORE INTO " + _s_tmp_region_tbl + " "
			"(label, chr, posMin, posMax) VALUES (?,?,?,?)";

	sqlite3_stmt* insert_stmt;
	sqlite3_prepare_v2(db, insert_sql.c_str(), -1, &insert_stmt, NULL);

	// Find the file and upload it...
	ifstream data_file(fn.c_str());
	if (!data_file.is_open()) {
		std::cerr << "WARNING: cannot find " << fn << ", ignoring.";
	} else {
		string line;
		vector<string> result;
		while (data_file.good()) {
			getline(data_file, line);
			split(result, line, is_any_of(" \n\t"), boost::token_compress_on);
			if(result.size() == 4 && result[0][0] == '#'){
				short chr = Locus::getChrom(result[0]);
				int posMin = atoi(result[2].c_str());
				int posMax = atoi(result[3].c_str());

				if(chr != -1 && posMin && posMax){
					sqlite3_bind_int(insert_stmt, 0, chr);
					sqlite3_bind_int(insert_stmt, 2, posMin);
					sqlite3_bind_int(insert_stmt, 3, posMax);
					sqlite3_bind_text(insert_stmt, 1, result[1].c_str(), -1, SQLITE_STATIC);

					while(sqlite3_step(insert_stmt)==SQLITE_ROW) {}

					sqlite3_reset(insert_stmt);
				}
			}
		}
		data_file.close();
	}
}

void RegionCollectionSQLite::updateZones(){
	// See loki_updater.py for the algorithm here

	// Get the zone size
	int zone_size=100000;
	string zone_sql = "SELECT value FROM setting "
			"WHERE setting='biopolymer_zone_size'";
	sqlite3_exec(db, zone_sql.c_str(), &parseSingleIntQuery, &zone_size, NULL);

	// Reverse any regions that are backwards
	string region_reverse_sql = "UPDATE " + _s_tmp_region_tbl + " "
			"SET posMin = posMax, posMax = posMin WHERE posMin > posMax";
	sqlite3_exec(db, region_reverse_sql.c_str(), NULL, NULL, NULL);

	// Find the minimum and maximum zone sizes
	int minPos, maxPos = 0;
	string min_pos_sql = "SELECT MIN(posMin) FROM " + _s_tmp_region_tbl;
	string max_pos_sql = "SELECT MAX(posMax) FROM " + _s_tmp_region_tbl;

	sqlite3_exec(db, min_pos_sql.c_str(), &parseSingleIntQuery, &minPos, NULL);
	sqlite3_exec(db, max_pos_sql.c_str(), &parseSingleIntQuery, &maxPos, NULL);

	minPos = minPos / zone_size;
	maxPos = maxPos / zone_size;

	// Insert all zones needed into a temporary table
	string tmp_zone_tbl = "__t_zone_tmp";
	string zone_tbl_sql = "CREATE TEMPORARY TABLE " + tmp_zone_tbl + " "
			"(zone INTEGER PRIMARY KEY NOT NULL)";
	sqlite3_exec(db, zone_tbl_sql.c_str(), NULL, NULL, NULL);

	string zone_insert_query = "INSERT INTO " + tmp_zone_tbl + " VALUES (?)";
	sqlite3_stmt* zone_ins_stmt;
	sqlite3_prepare_v2(db, zone_insert_query.c_str(), -1, &zone_ins_stmt, NULL);

	for (int i=minPos; i <= maxPos; i++){
		sqlite3_bind_int(zone_ins_stmt, 1, i);
		while(sqlite3_step(zone_ins_stmt) == SQLITE_ROW){}
		sqlite3_reset(zone_ins_stmt);
	}
	sqlite3_finalize(zone_ins_stmt);

	// OK, now we're ready to do some insertin'!
	stringstream zone_ins_str;
	zone_ins_str << "INSERT OR IGNORE INTO " << _s_tmp_zone_tbl << " (region_id,chr,zone) "
			<< "SELECT rb.biopolymer_id, rb.chr, z.zone "
			<< "FROM " << _s_tmp_region_tbl << " AS rb JOIN " << tmp_zone_tbl << " AS z "
			<< "ON z.zone >= rb.posMin / " << zone_size << " "
			<< "AND z.zone <= rb.posMax / " << zone_size << " ";

	string zone_ins_sql = zone_ins_str.str();
	sqlite3_exec(db, zone_ins_sql.c_str(), NULL, NULL, NULL);

	string tmp_tbl_drop = "DROP TABLE " + tmp_zone_tbl;
	sqlite3_exec(db, tmp_tbl_drop.c_str(), NULL, NULL, NULL);

	string index_sql = "CREATE INDEX " + _s_tmp_zone_tbl + "__zone ON " +
			_s_tmp_zone_tbl + " (chr, zone, region_id)";
	sqlite3_exec(db, index_sql.c_str(), NULL, NULL, NULL);
}

uint RegionCollectionSQLite::Load(const unordered_set<uint>& ids,
		const vector<string>& alias_list){

	loadFiles();

	// First things first, get a list of all of the ids associated with aliases
	unordered_set<uint> id_list(ids);
	if (alias_list.size() > 0){
		stringstream alias_stream;
		vector<string>::const_iterator a_itr = alias_list.begin();

		alias_stream << "SELECT biopolymer_id FROM biopolymer_name WHERE biopolymer_name IN ('"
				<< *a_itr << "'";

		while(++a_itr != alias_list.end()){
			alias_stream << ",'" << *a_itr << "'";
		}

		alias_stream << ")";

		sqlite3_exec(db, alias_stream.str().c_str(), parseRegionIDQuery, &id_list, NULL);
	}

	string where_clause = "WHERE ";

	if (id_list.size() > 0) {
		stringstream id_stream;
		unordered_set<uint>::const_iterator itr = id_list.begin();

		id_stream << "biopolymer.biopolymer_id IN (" << *itr;
		while(++itr != id_list.end()){
			id_stream << "," << *itr;
		}
		id_stream << ") AND ";
		where_clause += id_stream.str();
	}

	string command = "SELECT biopolymer.biopolymer_id, biopolymer.label "
			"FROM biopolymer_zone "
			"INNER JOIN biopolymer_region USING (biopolymer_id, chr) "
			"INNER JOIN biopolymer USING (biopolymer_id) ";

	where_clause += "biopolymer_zone.chr=:chrom AND biopolymer_region.ldprofile_id IN (:pop_id, :def_pop_id) "
			"AND zone=:pos_zone AND biopolymer_region.posMin<=:pos AND biopolymer_region.posMax>=:pos ";

	string stmt = command + where_clause;

	sqlite3_stmt* region_stmt;

	sqlite3_prepare_v2(db, stmt.c_str(), -1, &region_stmt, NULL);

	int pop_idx = sqlite3_bind_parameter_index(region_stmt, ":pop_id");
	int def_pop_idx = sqlite3_bind_parameter_index(region_stmt, ":def_pop_id");
	int chr_idx = sqlite3_bind_parameter_index(region_stmt, ":chrom");
	int pos_zone_idx = sqlite3_bind_parameter_index(region_stmt, ":pos_zone");
	int pos_idx = sqlite3_bind_parameter_index(region_stmt, ":pos");

	sqlite3_bind_int(region_stmt, pop_idx, _popID);
	sqlite3_bind_int(region_stmt, def_pop_idx, _def_id);

	sqlite3_stmt* tmp_region_stmt;

	string tmp_region_sql = "SELECT region_id, label, chr, posMin, posMax "
			"FROM " + _s_tmp_zone_tbl + " AS z "
			"INNER JOIN " + _s_tmp_region_tbl + " AS r USING (region_id, chr) "
			"WHERE zone=:pos_zone AND z.chr=:chrom AND posMin<=:pos AND posMax>=:pos";

	sqlite3_prepare_v2(db, tmp_region_sql.c_str(), -1, &tmp_region_stmt, NULL);

	int tmp_pos_zone_idx = sqlite3_bind_parameter_index(tmp_region_stmt, ":pos_zone");
	int tmp_chr_idx = sqlite3_bind_parameter_index(tmp_region_stmt, ":chrom");
	int tmp_pos_idx = sqlite3_bind_parameter_index(tmp_region_stmt, ":pos");


	// OK, now iterate over the dataset we have
	Container::const_iterator itr = _dataset->begin();
	int zone_size = _info->getZoneSize();
	while(itr != _dataset->end()){
		sqlite3_bind_int(region_stmt, chr_idx, (*itr)->getChrom());
		sqlite3_bind_int(region_stmt, pos_idx, static_cast<int>((*itr)->getPos()));
		sqlite3_bind_int(region_stmt, pos_zone_idx, (*itr)->getPos()/zone_size);

		sqlite3_bind_int(tmp_region_stmt, tmp_chr_idx, (*itr)->getChrom());
		sqlite3_bind_int(tmp_region_stmt, tmp_pos_idx, static_cast<int>((*itr)->getPos()));
		sqlite3_bind_int(tmp_region_stmt, tmp_pos_zone_idx, (*itr)->getPos()/zone_size);

		// Now, execute the query!
		while(sqlite3_step(region_stmt) == SQLITE_ROW){
			Knowledge::Region* reg = addRegion(region_stmt);
			reg->addLocus(*itr);
			_locus_map[*itr].insert(reg);
		}
		sqlite3_reset(region_stmt);

		while(sqlite3_step(tmp_region_stmt) == SQLITE_ROW){
			uint pop_id_result = sqlite3_column_int(tmp_region_stmt, 0);
			string label = (const char*) (sqlite3_column_text(tmp_region_stmt, 1));
			short chr = static_cast<short>(sqlite3_column_int(tmp_region_stmt, 2));
			uint start = static_cast<uint>(sqlite3_column_int(tmp_region_stmt, 3));
			uint end = static_cast<uint>(sqlite3_column_int(tmp_region_stmt, 4));
			Knowledge::Region* reg = AddRegion(label, pop_id_result, chr, start, end);
			reg->addLocus(*itr);
			_locus_map[*itr].insert(reg);
		}
		sqlite3_reset(tmp_region_stmt);

		++itr;
	}

	sqlite3_finalize(region_stmt);

	return _region_map.size();

}

void RegionCollectionSQLite::prepareStmts(){

	string region_alias_sql = "SELECT name "
			"FROM biopolymer_name WHERE biopolymer_id=?";

	sqlite3_prepare_v2(db, region_alias_sql.c_str(), -1, &_region_name_stmt, NULL);

	string region_bound_sql = "SELECT ldprofile_id, chr, posMin, posMax "
			"FROM biopolymer_bound "
			"WHERE biopolymer_id=? AND ldprofile_id IN (?, ?)";

	sqlite3_prepare_v2(db, region_bound_sql.c_str(), -1, &_region_bound_stmt, NULL);

	_popID = _info->getPopulationID(pop_str);
	_def_id = _info->getPopulationID("n/a");

	sqlite3_bind_int(_region_bound_stmt, 2, _popID);
	sqlite3_bind_int(_region_bound_stmt, 3, _def_id);

}

int RegionCollectionSQLite::parseRegionIDQuery(void* obj, int ncols, char** colVals, char** colNames){
	unordered_set<uint>* id_list = static_cast<unordered_set<uint>* >(obj);

	if (ncols > 1){
		// TOo many columns!
		return 2;
	}

	uint id=atoi(colVals[0]);
	if(id){
		id_list->insert(id);
	}
	return 0;

}

Knowledge::Region* RegionCollectionSQLite::addRegion(sqlite3_stmt* row){
	uint id = static_cast<uint>(sqlite3_column_int(row, 0));
	if(_region_map.find(id) != _region_map.end()){
		return (* _region_map.find(id)).second;
	} else {
		string name = (const char*) sqlite3_column_text(row, 1);
		Region* new_reg = new Region(name, id);

		sqlite3_bind_int(_region_name_stmt, 1, id);

		while(sqlite3_step(_region_name_stmt) == SQLITE_ROW){
			new_reg->addAlias((const char *) sqlite3_column_text(_region_name_stmt, 0));
		}
		sqlite3_reset(_region_name_stmt);

		sqlite3_bind_int(_region_bound_stmt, 1, id);
		while(sqlite3_step(_region_bound_stmt) == SQLITE_ROW){
			int pop_id_result = sqlite3_column_int(_region_bound_stmt, 0);
			short chr = static_cast<short>(sqlite3_column_int(_region_bound_stmt, 1));
			uint start = static_cast<uint>(sqlite3_column_int(_region_bound_stmt, 2));
			uint end = static_cast<uint>(sqlite3_column_int(_region_bound_stmt, 3));

			if(pop_id_result == _def_id){
				new_reg->addDefaultBoundary(chr, start, end);
			}else if(pop_id_result == _popID){
				new_reg->addPopulationBoundary(chr, start, end);
			}
		}
		sqlite3_reset(_region_bound_stmt);

		insertRegion(*new_reg);

		return new_reg;
	}
}

int RegionCollectionSQLite::parseSingleIntQuery(void* pop_id, int n_cols, char** col_vals, char** col_names){
	if(n_cols !=  1){
		return 2;
	}

	int* result = (int*) pop_id;
	(*result) = atoi(col_vals[0]);
	return 0;
}

} // namespace Knowledge;
