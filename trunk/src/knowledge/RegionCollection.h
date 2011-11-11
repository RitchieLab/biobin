#ifndef KNOWLEDGE_REGIONCOLLECTION_H
#define KNOWLEDGE_REGIONCOLLECTION_H

#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "knowledge/region.h"

using boost::unordered_map;
using std::string;
using std::vector;
using boost::unordered_set;

namespace Knowledge{

/**
 * This class is a complete rewrite of the RegionManager class in 
 * regionmanager.h.  It is designed to eliminate the id->index confusion.  
 */
class RegionCollection{
	
public:
	RegionCollection():region_not_found("Not Found",-1,-1){};
	virtual ~RegionCollection(){};

	/**
	 * Adds region and returns
    * @param name		the primary name
    * @param id		the gene_id from the database
    * @param start	assigned to both true and effective start
    * @param stop		assigned to both true and effective stop
	 * @param aliases Comma separated list of aliases
    * @return returns a reference to the actual region
    */
	void AddRegion(const char *, uint id, char chrom, uint start, uint stop, const char *aliases = "");

	/**
	 * Adds region and returns it's index
    * @param name		primary name
    * @param id		gene_id
    * @param effStart	effective Start
    * @param effStop		effective stop
    * @param trueStart	true start
    * @param trueStop	true stop
	 * @param aliases Comma separated list of aliases
    * @return returns a reference to the actual region
    */
	void AddRegion(const char *name, uint id, char chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const char* aliases = "");
	
	/**
	 * Access bracket operator, indexing by id
	 * NOTE: contrary to the default behavior of a map, this will not insert
	 * values into the map.  To insert, you must use AddRegion
	 */
	Region& operator[](const uint);

	/**
	 * Access bracket operator, indexing by alias
	 * NOTE: contrary to the default behavior of a map, this will not insert
	 * values into the map.  To insert, you must use AddRegion
	 */
	Region& operator[](const string&);
	/**
	 * Access bracket operator, indexing by ID.
	 */	
	const Region& operator[](const uint) const;
	/** 
	 * Access bracket operator, indexing by name.  NOTE: there is no 
	 * assignment bracket operator by name/alias. 
	 */
	const Region& operator[](const string&) const;
	
	/**
	 * Loading function - must be subclassed
	 * This function does the heavy lifting, and is the only function that needs
	 * to be subclassed
	 */
	virtual uint Load(const uint popID,
			const unordered_set<uint>& ids,
			const vector<string>& aliasList) = 0;

	/**
	 * Calls Load(...) with an empty ID list
	 */
	virtual uint Load(const uint popID,
			const vector<string>& aliasList);

	/**
	 * Calls Load(...) with an empty aliasList
	 */
	virtual uint Load(const uint popID,
			const unordered_set<uint>& ids);

	/**
	 * Calls Load(...) with an empty set of ids and and empty aliasList
	 */
	virtual uint Load(const uint popID);

	/**
	 * Calls Load(...) with popID=0, empty set of ids and empty aliasList
	 */
	virtual uint Load();

	/**
	 * Compares the region agains the private flagged region
	 */
	bool isValid(const Region&);

	/**
	 * Special value used for if the requested region is not in the map
	 */

protected:
	// A map from id -> Region
	unordered_map<uint,Region> region_map;
	// A map from alias -> id (to get region, essentially use
	// regionMap[aliasMap[alias]])
	unordered_map<string,uint> alias_map;

	/**
	 * Deletes regions that have no SNPs associated with them
	 */
	void Squeeze();

private:
	/**
	 * Copy constructor + assignment operator.  This object should not be
	 * copied or assinged.  EVER!
	 */
	RegionCollection(const RegionCollection&);
	RegionCollection& operator=(const RegionCollection&);

	Region region_not_found;

};
}


#endif
