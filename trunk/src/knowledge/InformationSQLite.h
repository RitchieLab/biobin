/*
 * InformationSQLite.h
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_INFORMATIONSQLITE_H
#define KNOWLEDGE_INFORMATIONSQLITE_H

#include <string>

#include <sqlite3.h>

#include "Information.h"

using std::string;

namespace Knowledge{

class InformationSQLite : public Information {
public:
	virtual ~InformationSQLite();

	InformationSQLite(const string& filename);
	InformationSQLite(sqlite3* db);

	virtual int getPopulationID(const string& pop_str);
	virtual const string getResourceVersion(const string& resource);

	virtual void getGroupTypes(const vector<int>& type_ids,
			map<int, string>& group_types_out);

private:
	// SQLite callback to parse a single column, which is returned via the
	// first argument, which must be a pointer to a string.
	static int parseSingleStringQuery(void*, int, char**, char**);

	static int parseGroupTypeQuery(void*, int, char**, char**);

	sqlite3* _db;
	bool _self_open;

};

}




#endif /* INFORMATIONSQLITE_H_ */
