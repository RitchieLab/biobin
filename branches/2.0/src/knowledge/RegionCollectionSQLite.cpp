/*
 * RegionCollectionSQLite.cpp
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#include "RegionCollectionSQLite.h"

// Use the straight-up sqlite interface
#include <sqlite3.h>
#include <stdlib.h>
#include <sstream>

#include "InformationSQLite.h"

#include "Locus.h"

using std::stringstream;

namespace Knowledge{

RegionCollectionSQLite::RegionCollectionSQLite(const string& fn) : self_open(true){
	sqlite3_open(fn.c_str(),&db);
	_info = new InformationSQLite(db);
}

RegionCollectionSQLite::RegionCollectionSQLite(sqlite3* db_conn) :
		self_open(false), db(db_conn) {
	_info = new InformationSQLite(db);
}

RegionCollectionSQLite::~RegionCollectionSQLite(){
	if (self_open){
		sqlite3_close(db);
	}
	delete _info;
}
uint RegionCollectionSQLite::Load(const unordered_set<uint>& ids,
		const vector<string>& alias_list){

	string where_clause = "";

	if (ids.size() > 0) {
		stringstream id_stream;
		unordered_set<uint>::const_iterator itr = ids.begin();
		unordered_set<uint>::const_iterator end = ids.end();
		id_stream << "regions.gene_id IN (" << *itr;
		while(++itr != end){
			id_stream << "," << *itr;
		}
		id_stream << ")";
		where_clause = string(" WHERE ") + id_stream.str();
	}
	if (alias_list.size() > 0){
		stringstream alias_stream;
		vector<string>::const_iterator itr = alias_list.begin();
		vector<string>::const_iterator end = alias_list.end();
		alias_stream << "region_aliases IN (" << *itr;
		while(++itr != end){
			alias_stream << "," << *itr;
		}
		alias_stream << ")";
		if (ids.size() > 0){
			where_clause += string("OR ") + alias_stream.str();
		}
		else{
			where_clause = string(" WHERE ") + alias_stream.str();
		}
	}

	int popID = _info->getPopulationID(pop_str);

	string bound_vals = "DefBounds.start, DefBounds.end";
	string bound_tables = "LEFT OUTER JOIN region_bounds DefBounds "
			"ON regions.gene_id=DefBounds.gene_id AND "
			"DefBounds.population_id=0";
	if (popID != 0){
		bound_vals += " PopBounds.start, PopBounds.end";
		bound_tables += " LEFT OUTER JOIN region_bounds PopBounds "
				"ON regions.gene_id=PopBounds.gene_id AND "
				"PopBounds.population_id=";
		bound_tables += popID;
	}
	string command = string("SELECT regions.gene_id, chrom, primary_name, ") +
			bound_vals + string(", group_concat(alias) FROM regions ") +
			bound_tables + string(" LEFT OUTER JOIN region_alias ON "
				"regions.gene_id=region_alias.gene_id") + where_clause +
				string(" GROUP BY regions.gene_id");

	return sqlite3_exec(db, command.c_str(), parseRegionQuery, this, NULL);
}

int RegionCollectionSQLite::parseRegionQuery(void* obj, int ncols, char** colVals, char** colNames){
	RegionCollection* regions = (RegionCollection*) obj;

	if(ncols < 6){
		// Not enough columns!!
		return 2;
	}

	int id_idx = 0;
	int chrom_idx = 1;
	int name_idx = 2;
	int start_idx = 3;
	int end_idx = 4;
	int alias_idx = 5;

	// The number of additional columns selected in a population-specific query
	int pop_offset = 2;

	int id = atoi(colVals[id_idx]);
	int start = atoi(colVals[start_idx]);
	int end = atoi(colVals[end_idx]);

	short chrom = Locus::getChrom(colVals[chrom_idx]);

	// In this case, no population ID was given
	if (ncols == 6){
		regions->AddRegion(colVals[name_idx], id, chrom, start, end,
				start-gene_expansion, end+gene_expansion, colVals[alias_idx]);
	}else if (ncols == 8){
		int pop_start = colVals[start_idx + pop_offset] ? atoi(colVals[start_idx + pop_offset]) : start;
		int pop_end = colVals[end_idx + pop_offset] ? atoi(colVals[end_idx + pop_offset]) : end;
		regions->AddRegion(colVals[name_idx], id, chrom, start, end,
				pop_start, pop_end, colVals[alias_idx + pop_offset]);
	}else{
		// WE SHOULD NOT BE HERE!!
		return 1;
	}
	return 0;
}

} // namespace Knowledge;
