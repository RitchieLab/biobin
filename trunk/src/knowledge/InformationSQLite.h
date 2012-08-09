/*
 * InformationSQLite.h
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_INFORMATIONSQLITE_H
#define KNOWLEDGE_INFORMATIONSQLITE_H

#include <string>
#include <stdlib.h>
#include <map>

#include <sqlite3.h>

#include "Information.h"

using std::string;
using std::map;

namespace Knowledge{

/*!
 * \brief SQLite backend for the Information class.
 * Implements the Information interface using SQLite backend of the LOKI
 * database.
 */
class InformationSQLite : public Information {
public:
	virtual ~InformationSQLite();

	/*!
	 * \brief Creates an Infomration object with a filename.
	 * Takes the filename of the SQLite database and creates a Information
	 * object.
	 *
	 * \param filename The SQLite file containing the LOKI database
	 */
	InformationSQLite(const string& filename);
	/*!
	 * \brief Creates an Information object with an sqlite DB connection.
	 * Creates an Information object using an already open connection to a
	 * LOKI database through an sqlite object.  The database connection must
	 * remain open throughout the lifetime of this object.
	 *
	 * \param db An open SQLite DB object to a LOKI database.
	 */
	InformationSQLite(sqlite3* db);

	virtual int getPopulationID(const string& pop_str);
	virtual void getGroupTypes(const set<uint>& type_ids,
			map<int, string>& group_types_out);
	virtual int getZoneSize();

	virtual int getSNPRole(const Locus& loc, const Region& reg);

	virtual void printPopulations(ostream& os);
	virtual void printSources(ostream& os);

	virtual void loadRoles();

	virtual const set<unsigned int>& getSourceIds();

private:
	void prepRoleStmt();
	void prepRoleTables();

	static string _role_region_tbl;
	static string _role_zone_tbl;

	// blatantly stolen from ldsplineimporter!
	void getAndDropIndexes(const string& tbl_name, map<string, string>& indexes_out);
	void restoreIndexes(const string& tbl_name, const map<string, string>& index_map);
	void UpdateZones();

	static int parseRegionIndex(void*, int, char**, char**);


	/*!
	 *  SQLite callback to parse a single column, which is returned via the
	 * first argument, which must be a pointer to a string.
	 */
	static int parseSingleStringQuery(void*, int, char**, char**);

	/*!
	 * SQLite callback to parse the group type.  The first argument must be
	 * a map of integers to strings.
	 */
	static int parseGroupTypeQuery(void*, int, char**, char**);

	/*!
	 * SQLite callback to parse the zone
	 */
	static int parseSingleIntQuery(void*, int, char**, char**);
	static int parseMultiIntQuery(void*, int, char**, char**);
	static int printQueryResult(void*, int, char**, char**);


	sqlite3* _db;
	bool _self_open;

	// A mapping of db integer roles to SNP roles
	map<int, Information::snp_role> _role_map;

	// The prepared statement for getting SNP roles
	sqlite3_stmt* _role_stmt;

	// Prepared statment for getting SNP roles from regions
	sqlite3_stmt* _region_role_stmt;

};

}




#endif /* INFORMATIONSQLITE_H_ */
