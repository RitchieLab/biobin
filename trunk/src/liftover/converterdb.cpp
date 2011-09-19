/* 
 * File:   converterdb.cpp
 * Author: torstees
 * 
 * Created on May 17, 2011, 10:34 AM
 */

#include "converterdb.h"

namespace LiftOver {

int ConverterDB::LoadFromDB(const char *orig, soci::session& sociDB) {
	std::string build;
	std::string origBuild(orig);
	sociDB<<"SELECT version FROM versions WHERE element='build'", soci::into(build);
	newBuild = build;

	// 0 Is legitimate if we are already at the right build
	uint count = 0;
	
	if (origBuild != newBuild) {
		soci::rowset<soci::row> rs = (sociDB.prepare<<"SELECT chain_data FROM chain_files NATURAL JOIN build_versions WHERE build = :build", soci::use(origBuild));
		for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			soci::row const& row = *itr;
			std::string chain		= row.get<std::string>(0);
			AddChain(chain.c_str());
			count += 1;
		}
		
		//If this is true, then we didn't find a valid build, which is an error
		if (count == 0)
			return -1;
	}
	return count;
}

}

