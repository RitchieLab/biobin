/*
 * RegionCollectionSQLite.cpp
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#include "RegionCollectionSQLite.h"

#include <stdlib.h>
#include <sstream>


#include "Locus.h"

using std::stringstream;

namespace Knowledge{


RegionCollectionSQLite::~RegionCollectionSQLite(){
	if (self_open){
		sqlite3_close(db);
	}
}

uint RegionCollectionSQLite::Load(const unordered_set<uint>& ids,
		const vector<string>& alias_list){

	// First things first, get a list of all of the ids associated with aliases
	unordered_set<uint> id_list(ids);
	if (alias_list.size() > 0){
		stringstream alias_stream;
		vector<string>::const_iterator a_itr = alias_list.begin();

		alias_stream << "SELECT region_id FROM region_name WHERE region_name IN ('"
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

		id_stream << "region.region_id IN (" << *itr;
		while(++itr != id_list.end()){
			id_stream << "," << *itr;
		}
		id_stream << ") AND ";
		where_clause += id_stream.str();
	}

	int popID = _info->getPopulationID(pop_str);
	int def_id = _info->getPopulationID("n/a");

	string command = "SELECT region.region_id, region_bound.chr, region.label, "
			"region_bound.posMin, region_bound.posMax, rb_default.posMin, "
			"rb_default.posMax, group_concat(region_name.name) "
			"FROM region_bound "
			"INNER JOIN region_zone USING(region_id, chr, population_id) "
			"INNER JOIN region USING (region_id) "
			"LEFT JOIN region_bound AS rb_default "
				"ON region_bound.region_id=rb_default.region_id "
					"AND rb_default.population_id=:def_pop_id "
			"LEFT JOIN region_name ON region.region_id=region_name.region_id ";

	where_clause += "region_bound.chr=:chrom AND region_bound.population_id=:pop_id "
			"AND zone=:pos_zone AND region_bound.posMin<:pos AND region_bound.posMax>:pos ";

	string stmt = command + where_clause + "GROUP BY region.region_id";

	sqlite3_stmt* region_stmt;

	sqlite3_prepare_v2(db, stmt.c_str(), -1, &region_stmt, NULL);

	int pop_idx = sqlite3_bind_parameter_index(region_stmt, ":pop_id");
	int def_pop_idx = sqlite3_bind_parameter_index(region_stmt, ":def_pop_id");
	int chr_idx = sqlite3_bind_parameter_index(region_stmt, ":chrom");
	int pos_zone_idx = sqlite3_bind_parameter_index(region_stmt, ":pos_zone");
	int pos_idx = sqlite3_bind_parameter_index(region_stmt, ":pos");

	sqlite3_bind_int(region_stmt, pop_idx, popID);
	sqlite3_bind_int(region_stmt, def_pop_idx, def_id);

	// OK, now iterate over the dataset we have
	Container::const_iterator itr = _dataset->begin();
	int zone_size = _info->getZoneSize();
	while(itr != _dataset->end()){
		sqlite3_bind_int(region_stmt, chr_idx, (*itr)->getChrom());
		sqlite3_bind_int(region_stmt, pos_idx, static_cast<int>((*itr)->getPos()));
		sqlite3_bind_int(region_stmt, pos_zone_idx, (*itr)->getPos()/zone_size);

		// Now, execute the query!
		while(sqlite3_step(region_stmt) == SQLITE_ROW){
			Knowledge::Region* reg = addRegion(region_stmt);
			reg->addLocus(*itr);
		}
		sqlite3_reset(region_stmt);
		++itr;
	}

	sqlite3_finalize(region_stmt);

	return _region_map.size();

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
		return (*_region_map.find(id)).second;
	} else {
		short chrom = static_cast<short>(sqlite3_column_int(row, 1));
		string name = (const char*) sqlite3_column_text(row, 2);
		uint eff_start = static_cast<uint>(sqlite3_column_int(row,3));
		uint eff_end = static_cast<uint>(sqlite3_column_int(row,4));
		uint def_start = static_cast<uint>(sqlite3_column_int(row,5));
		uint def_end = static_cast<uint>(sqlite3_column_int(row,6));
		string aliases = (const char*) sqlite3_column_text(row,7);

		return AddRegion(name, id, chrom, eff_start, eff_end, def_start, def_end, aliases);
	}
}

} // namespace Knowledge;
