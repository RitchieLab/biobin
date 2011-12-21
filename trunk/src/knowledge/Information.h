/*
 * Information.h
 *
 *  Created on: Dec 2, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_INFORMATION_H
#define KNOWLEDGE_INFORMATION_H

#include <string>
#include <ostream>
#include <vector>
#include <map>
#include <set>

using std::ostream;
using std::string;
using std::vector;
using std::map;
using std::set;

namespace Knowledge{

/**
 * Defines a class that gets general information about the database
 */
class Information{
public:

	virtual ~Information(){}

	virtual int getPopulationID(const string& pop_str) = 0;
	virtual const string getResourceVersion(const string& resource) = 0;

	// The following models are abstractions of junk that I found in Biofilter's
	// main application class.  See below for the biofilter code for these
	// functions
	virtual void listPopulationIDs(ostream& os){}
	virtual void listGroupIDs(ostream& os, const vector<string>& group_list){}
	virtual void listRegions(ostream& os, const vector<string>& alias_list, const vector<string>& alias_type){}

	virtual void getGroupTypes(const set<uint>& group_ids,
			map<int, string>& type_names_out) = 0;
};
/* listRegions
std::string aliasClause = "";
std::string typeClause = "";
if (aliasList.size() > 0)
	aliasClause = std::string("( alias LIKE '") + Utility::Join(aliasList, "%' OR alias LIKE '") + "%' ) ";
if (aliasType.size() > 0)
	typeClause = std::string("( region_alias_type_desc LIKE '") + Utility::Join(aliasType, "%' OR region_alias_type_desc LIKE '") + "%' )";
std::string whereClause = aliasClause + typeClause;
if (aliasClause.length() > 0 && typeClause.length() > 0)
	whereClause = std::string("AND ")  + aliasClause + " AND " + typeClause;
else if (whereClause.length() > 0)
	whereClause = std::string("AND ") + whereClause;
os<<"Name\tAlias\tChrom\tStart\tEnd\tDescription\tAlias Type\n";
soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT primary_name, alias, chrom, start, end, description, region_alias_type_desc FROM regions NATURAL JOIN region_bounds NATURAL JOIN region_alias NATURAL JOIN region_alias_type WHERE population_id=0 "<<whereClause);
for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
	soci::row const& row = *itr;
	os<<row.get<std::string>(0)<<"\t"
		 <<row.get<std::string>(1)<<"\t"
		 <<row.get<std::string>(2)<<"\t"
		 <<row.get<int>(3)<<"\t"
		 <<row.get<int>(4)<<"\t"
		 <<row.get<std::string>(5)<<"\t"
		 <<row.get<std::string>(6)<<"\n";
}*/

/* listGroupIDs
std::string clause = "";

if (searchList.size() > 0)
	clause = std::string("WHERE ( group_name LIKE '%") + Utility::Join(searchList, "%' OR group_name LIKE '%") + "%' OR group_desc LIKE '%" + Utility::Join(searchList, "%' OR group_desc LIKE '%'") + "%') ";
os<<"ID\tName\tDescription\n";
soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT group_id, group_name, group_desc FROM groups "<<clause);
for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
	soci::row const& row = *itr;
	os<<row.get<int>(0)<<"\t"
		 <<row.get<std::string>(1)<<"\t"
		 <<row.get<std::string>(2)<<"\n";
}
*/
/* listPopulationIDs
soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT population_label, pop_description FROM populations");
os<<"Label\tDescription\n";
for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
	soci::row const& row = *itr;
	os<<row.get<std::string>(0)<<"\t"<<row.get<std::string>(1)<<"\n";
}
*/
}


#endif /* INFORMATION_H_ */
