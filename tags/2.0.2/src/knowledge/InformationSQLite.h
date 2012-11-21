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
	InformationSQLite(const std::string& filename);
	/*!
	 * \brief Creates an Information object with an sqlite DB connection.
	 * Creates an Information object using an already open connection to a
	 * LOKI database through an sqlite object.  The database connection must
	 * remain open throughout the lifetime of this object.
	 *
	 * \param db An open SQLite DB object to a LOKI database.
	 */
	InformationSQLite(sqlite3* db);

	virtual int getPopulationID(const std::string& pop_str);
	virtual void getGroupTypes(const std::set<unsigned int>& type_ids,
			std::map<int, std::string>& group_types_out);
	virtual int getZoneSize();

	virtual unsigned long getSNPRole(const Locus& loc, const Region& reg, bool use_cache);

	virtual void printPopulations(std::ostream& os);
	virtual void printSources(std::ostream& os);

	virtual void loadRoles(const RegionCollection& reg);

	virtual const std::set<unsigned int>& getSourceIds();

private:
	void prepRoleStmt();
	void prepRoleTables();

	static std::string _role_region_tbl;
	static std::string _role_zone_tbl;

	// blatantly stolen from ldsplineimporter!
	void getAndDropIndexes(const std::string& tbl_name, std::map<std::string, std::string>& indexes_out);
	void restoreIndexes(const std::string& tbl_name, const std::map<std::string, std::string>& index_map);
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
	std::map<int, Information::snp_role> _role_map;

	// The prepared statement for getting SNP roles
	sqlite3_stmt* _role_stmt;

	// Prepared statment for getting SNP roles from regions
	sqlite3_stmt* _region_role_stmt;

	// Zone size used in the temp tables
	int _tmp_zone_sz;

};

}




#endif /* INFORMATIONSQLITE_H_ */
