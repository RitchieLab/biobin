/* 
 * File:   converterdb.h
 * Author: torstees
 *
 * Created on May 17, 2011, 10:34 AM
 */

#ifndef KNOWLEDGE_CONVERTERSQLITE_H
#define	KNOWLEDGE_CONVERTERSQLITE_H

#include "Converter.h"

#include <string>
#include <sqlite3.h>

namespace Knowledge {
namespace Liftover {

class ConverterSQLite : public Converter {
public:
	ConverterSQLite(const std::string& orig_build, const std::string& db_fn);
	ConverterSQLite(const std::string& orig_build, sqlite3* db);
	virtual ~ConverterSQLite();
	
	virtual int Load();

private:

	static int parseChains(void*, int, char**, char**);
	static int parseChainData(void*, int, char**, char**);

	sqlite3* _db;
	bool _self_open;

};

} // namespace Liftover
} // namespace Knowledge
#endif	/* CONVERTERDB_H */

