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

/*!
 * \brief Defines a class that gets general information about the database.
 * This class is an interface to get specific information about the database.
 * The information returned by this class is meta-information about the
 * database, like the genomic build and the list of populations available.
 */
class Information{
public:

	/*!
	 * Destroys the Information object
	 */
	virtual ~Information(){}

	/*!
	 * \brief Converts a population string into a Population ID.
	 * Takes the given pop_str and returns an integer ID needed for DB queries.
	 *
	 * \param pop_str The Population label.
	 *
	 * \return The integer index of the given population, or 0 if not found.
	 */
	virtual int getPopulationID(const string& pop_str) = 0;

	/*!
	 * \brief Gets the names of the sources for the given IDs
	 * Queries the database and returns a map of IDs to source name.
	 *
	 * \param group_ids A set of IDs to obtain sourse information for.  If empty,
	 * return source information for all available sources
	 * \param[out] type_names_out A mapping of ID->source name
	 */
	virtual void getGroupTypes(const set<uint>& group_ids,
			map<int, string>& type_names_out) = 0;

	/*!
	 * \brief Gets the zone size.
	 * Queries the database and returns the size of the zones used in building
	 * the region_zone table.
	 */
	virtual int getZoneSize() = 0;

};

}


#endif /* INFORMATION_H_ */
