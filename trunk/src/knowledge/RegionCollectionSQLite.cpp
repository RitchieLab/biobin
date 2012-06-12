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
	sqlite3_finalize(_region_name_stmt);
	sqlite3_finalize(_region_bound_stmt);

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

	string command = "SELECT region.region_id, region.label "
			"FROM region_zone "
			"INNER JOIN region_bound USING (region_id, chr, population_id) "
			"INNER JOIN region USING (region_id) ";

	where_clause += "region_zone.chr=:chrom AND region_zone.population_id IN (:pop_id, :def_pop_id) "
			"AND zone=:pos_zone AND region_bound.posMin<:pos AND region_bound.posMax>:pos ";

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
			_locus_map[*itr].insert(reg);
		}
		sqlite3_reset(region_stmt);
		++itr;
	}

	sqlite3_finalize(region_stmt);

	return _region_map.size();

}

void RegionCollectionSQLite::prepareStmts(){

	string region_alias_sql = "SELECT name "
			"FROM region_name WHERE region_id=?";

	sqlite3_prepare_v2(db, region_alias_sql.c_str(), -1, &_region_name_stmt, NULL);

	string region_bound_sql = "SELECT population_id, chr, posMin, posMax "
			"FROM region_bound "
			"WHERE region_id=? AND population_id IN (?, ?)";

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

			//insertRegionBound(*new_reg, chr, start, end);

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

} // namespace Knowledge;
