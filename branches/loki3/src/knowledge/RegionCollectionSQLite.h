/*
 * RegionCollectionSQLite.h
 *
 *  Created on: Nov 10, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_REGIONCOLLECTIONSQLITE_H
#define KNOWLEDGE_REGIONCOLLECTIONSQLITE_H

#include "RegionCollection.h"
#include "InformationSQLite.h"

#include <sqlite3.h>

namespace Knowledge{

/*!
 * \brief A class that implements the region collection from a SQLite database.
 * This class represents a RegionCollection from an SQLite LOKI database.
 * Users may construct this from either an already open SQLite3 database
 * connection, or from a filename that gives the location of the SQLite3 db.
 * the only dependency required for this class is the sqlite3 C interface
 * library (v. >=3.5.6), which should come standard with most modern Linux
 * distributions.
 */
class RegionCollectionSQLite : public RegionCollection{

public:
	/*!
	 * \brief Create a RegionCollection giving the DB location.
	 * Creates a RegionCollection from the fiven location of the sqlite3
	 * database.  This file must exist and conform to the LOKI specifications.
	 *
	 * \param fn The filename of the LOKI database.
	 */
	template <class T_cont>
	RegionCollectionSQLite(const std::string& fn, const T_cont& loci, Information* info=0);
	/*!
	 * \brief Create a RegionCollection giving an open DB connection.
	 * Creates a RegionCollection with an already open databse connection that
	 * points to a LOKI database.  In this case, the user is responsible for
	 * closing the connection, and the connection cannot be closed until this
	 * object is destroyed.
	 *
	 * \param db The database connection to the LOKI database.
	 */
	template <class T_cont>
	RegionCollectionSQLite(sqlite3 *db, const T_cont& loci, Information* info=0);
	/*!
	 * Destroys the RegionCollection objects, and if necessary closes the
	 * database connection.
	 */
	virtual ~RegionCollectionSQLite();

	/*!
	 * \brief Loading function implementation.
	 * This function loads all of the Regions that are identified by the ids or
	 * by the aliasList from the LOKI database into memory.  If both ids and
	 * aliasList are empty, will load all Regions from the databse into memory.
	 *
	 * \param ids A list of IDs to filter the results by.
	 * \param aliasList A list of aliases to filter results
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	virtual unsigned int Load(const boost::unordered_set<unsigned int>& ids,
			const std::vector<std::string>& aliasList);

	virtual void loadFiles();

private:
	//! true if we opened the connection, false otherwise
	bool self_open;
	bool self_info;
	//! sqlite connection
	sqlite3 *db;

	sqlite3_stmt* _region_name_stmt;
	sqlite3_stmt* _region_bound_stmt;

	int _popID;
	int _def_id;
	int _schema_vers;

	static std::string _s_tmp_region_tbl;
	static std::string _s_tmp_zone_tbl;
	static std::string _s_tmp_name_tbl;
	static std::string _s_tmp_bound_tbl;

	//! Adds a region based on the row (or returns the already added region)
	Knowledge::Region* addRegion(sqlite3_stmt* row);

	//! Adds the stuff from a single region file
	void loadFile(const std::string& fn);

	//! Creates zones for the temp region table
	void updateZones();

	//! Prepares some statments to be used in finding the row information
	void prepareStmts();

	/*!
	 * Callback to parse a list of region IDs and add them to the undordered_list
	 */
	static int parseRegionIDQuery(void*, int, char**, char**);
	static int parseSingleIntQuery(void*, int, char**, char**);

};

template<class T_cont>
RegionCollectionSQLite::RegionCollectionSQLite(const std::string& fn,
		const T_cont& loci, Information* info) : RegionCollection(loci), self_open(true), self_info(!info) {
	sqlite3_open(fn.c_str(), &db);
	if(!info){
		_info = new InformationSQLite(db);
	}else{
		_info = info;
	}

	prepareStmts();
}

template<class T_cont>
RegionCollectionSQLite::RegionCollectionSQLite(sqlite3* db_conn,
		const T_cont& loci, Information* info) : RegionCollection(loci),
		self_open(false), self_info(!info), db(db_conn) {
	if(!info){
		_info = new InformationSQLite(db);
	}else{
		_info = info;
	}

	prepareStmts();
}

}




#endif /* REGIONCOLLECTIONSQLITE_H_ */
