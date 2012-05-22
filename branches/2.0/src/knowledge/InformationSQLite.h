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

#include <sqlite3.h>

#include "Information.h"

using std::string;

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

private:
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

	sqlite3* _db;
	bool _self_open;

};

}




#endif /* INFORMATIONSQLITE_H_ */
