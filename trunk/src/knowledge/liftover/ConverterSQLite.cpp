/* 
 * File:   converterdb.cpp
 * Author: torstees
 * 
 * Created on May 17, 2011, 10:34 AM
 */

#include "ConverterSQLite.h"

#include <sstream>
#include <iostream>

#include "Chain.h"


using std::stringstream;
using std::map;
using std::set;
using std::string;

namespace Knowledge {
namespace Liftover{

ConverterSQLite::ConverterSQLite(const string& orig_build,
		const string& db_filename) :
		Converter(orig_build), _self_open(true), _build_loaded(_origBuild + "_"){
	sqlite3_open(db_filename.c_str(),&_db);
}

ConverterSQLite::ConverterSQLite(const string& orig_build,
		sqlite3* db) :
		Converter(orig_build), _db(db), _self_open(false),  _build_loaded(_origBuild + "_"){}

ConverterSQLite::~ConverterSQLite(){
	if (_self_open){
		sqlite3_close(_db);
	}
}

int ConverterSQLite::Load() {

	if (_build_loaded != _origBuild) {
		_build_loaded = _origBuild;
		_chains.clear();
		// Find the current version that we are building to
		stringstream ss;
		ss << "SELECT chain_id, score, old_chr, old_start, old_end, new_chr, is_fwd "
				"FROM chain INNER JOIN grch_ucschg "
				"ON grch_ucschg.ucschg=chain.old_ucschg "
				"WHERE grch = '" << _origBuild << "' AND " <<
				"new_ucschg = (SELECT value FROM setting WHERE setting='ucschg');";

		sqlite3_exec(_db, ss.str().c_str(), parseChains, this, NULL);

		map<short, set<Chain*> >::iterator itr = _chains.begin();
		while(itr != _chains.end()){
			set<Chain*>::iterator s_itr = (*itr).second.begin();
			while(s_itr != (*itr).second.end()){
				stringstream ss_data;
				ss_data << "SELECT old_start, old_end, new_start "
						"FROM chain_data WHERE chain_id=" << (*s_itr)->getID() << ";";

				sqlite3_exec(_db, ss_data.str().c_str(), parseChainData, *s_itr, NULL);

				++s_itr;
			}

			++itr;
		}
	}

	return _chains.size();
}



int ConverterSQLite::parseChains(
		void* obj, int n_cols, char** col_vals, char** col_names){

	ConverterSQLite *conv = static_cast<ConverterSQLite*>(obj);

	// Not the right # of columns
	if (n_cols != 7){
		return 2;
	}

	int id = atoi(col_vals[0]);
	long score = atol(col_vals[1]);
	short o_chr = static_cast<short>(atoi(col_vals[2]));
	int o_start = atoi(col_vals[3]);
	int o_end = atoi(col_vals[4]);
	short n_chr = static_cast<short>(atoi(col_vals[5]));
	bool fwd = static_cast<bool>(atoi(col_vals[6]));

	conv->_chains[o_chr].insert(new Chain(id, score, o_start, o_end, n_chr, fwd));

	return 0;
}

int ConverterSQLite::parseChainData(
		void* obj, int n_cols, char** col_vals, char** col_names){

	Chain *chn = static_cast<Chain*>(obj);

	// Not the right number of columns
	if (n_cols != 3){
		return 2;
	}

	int old_s = atoi(col_vals[0]);
	int old_e = atoi(col_vals[1]);
	int new_s = atoi(col_vals[2]);

	chn->addSegment(old_s, old_e, new_s);

	return 0;
}

} // namespace Liftover
} // namespace Knowledge

