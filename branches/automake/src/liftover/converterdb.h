/* 
 * File:   converterdb.h
 * Author: torstees
 *
 * Created on May 17, 2011, 10:34 AM
 */

#ifndef CONVERTERDB_H
#define	CONVERTERDB_H

#include "converter.h"
#include <soci.h>
#include <string>
#include "utility/types.h"

namespace LiftOver {

class ConverterDB : public Converter {
public:
	ConverterDB();
	ConverterDB(const ConverterDB& orig);
	virtual ~ConverterDB();
	
	int LoadFromDB(const char *origBuild, soci::session& sociDB);
private:

};

inline
ConverterDB::ConverterDB() : Converter("", "") {}

inline
ConverterDB::ConverterDB(const ConverterDB& orig) : Converter(orig) {}

inline
ConverterDB::~ConverterDB() {}


}
#endif	/* CONVERTERDB_H */

