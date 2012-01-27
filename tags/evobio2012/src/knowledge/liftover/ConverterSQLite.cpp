/* 
 * File:   converterdb.cpp
 * Author: torstees
 * 
 * Created on May 17, 2011, 10:34 AM
 */

#include "ConverterSQLite.h"
#include <sqlite3.h>
#include <sstream>

using std::stringstream;

namespace Knowledge {
namespace Liftover{

ConverterSQLite::ConverterSQLite(const string& orig_build,
		const string& db_filename) :
		Converter(orig_build, orig_build), _self_open(true){
	sqlite3_open(db_filename.c_str(),&_db);
}

ConverterSQLite::ConverterSQLite(const string& orig_build,
		sqlite3* db) :
		Converter(orig_build, orig_build), _db(db), _self_open(false) {}

ConverterSQLite::~ConverterSQLite(){
	if (_self_open){
		sqlite3_close(_db);
	}
}

int ConverterSQLite::Load() {

	// Find the current version that we are building to

	int status =  sqlite3_exec(_db, "SELECT version FROM versions "
			"WHERE element='build'", parseCurrentVersion, this, NULL);

	if (status == 0){
		stringstream ss;
		ss << "SELECT chain_data FROM chain_files NATURAL JOIN build_versions "
				"WHERE build = '" << _newBuild << "';";

		status = sqlite3_exec(_db, ss.str().c_str(), parseChainFiles, this, NULL);
	}

	return status;
}

int ConverterSQLite::parseCurrentVersion(
		void* obj, int n_cols, char** col_vals, char** col_names){

	ConverterSQLite *conv = (ConverterSQLite*) obj;

	// Not the right number of columns
	if (n_cols != 1){
		return 2;
	}

	conv->_newBuild = col_vals[0];
	return 0;
}

int ConverterSQLite::parseChainFiles(
		void* obj, int n_cols, char** col_vals, char** col_names){

	ConverterSQLite *conv = (ConverterSQLite*) obj;

	// Not the right # of columns
	if (n_cols != 1){
		return 2;
	}

	conv->addChain(col_vals[0]);
	return 0;
}

} // namespace Liftover
} // namespace Knowledge

