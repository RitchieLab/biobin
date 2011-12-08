#ifndef KNOWLEDGE_REGIONCOLLECTION_H
#define KNOWLEDGE_REGIONCOLLECTION_H

#include <string>
#include <vector>
#include <set>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/icl/interval_map.hpp>

#include "Region.h"

using boost::unordered_map;
using boost::unordered_set;
using boost::icl::interval_map;
using boost::icl::interval;
using std::string;
using std::vector;
using std::set;


namespace Knowledge{

/**
 * This class is a complete rewrite of the RegionManager class in 
 * regionmanager.h.  It is designed to eliminate the id->index confusion.  
 */
class RegionCollection{
	
public:
	typedef set<Region*>::const_iterator const_region_iterator;

	RegionCollection():region_not_found("Not Found",-1,-1, 0, 0){};
	virtual ~RegionCollection();

	/**
	 * Adds region and returns
    * @param name		the primary name
    * @param id		the gene_id from the database
    * @param start	assigned to both true and effective start
    * @param stop		assigned to both true and effective stop
	 * @param aliases Comma separated list of aliases
    * @return returns a reference to the actual region
    */
	void AddRegion(const string& name, uint id, short chrom, uint start, uint stop, const string& aliases = "");

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
	void AddRegion(const string& name, uint id, short chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const string& aliases = "");
	
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
	//const set<Region&>& operator[](const string&);

	/**
	 * Access bracket operator, indexing by ID.
	 */	
	const Region& operator[](const uint) const;

	/** 
	 * Access by alias is controlled by an iterator
	 */
	const_region_iterator aliasBegin(const string&) const;
	const_region_iterator aliasEnd(const string&) const;


	const_region_iterator positionBegin(short chrom, uint pos) const;
	const_region_iterator positionEnd(short chrom, uint pos) const;
	/**
	 * Determine if a given alias is valid
	 */
	bool isValid(const string& alias) const {
		return _alias_map.find(alias) != _alias_map.end();
	}

	/**
	 * Compares the region agains the private flagged region
	 */
	bool isValid(const Region&) const;

	/**
	 * Takes a list of SNPs (realistically, just chrom->bp locations) and adds
	 * them to the Regions that contain them
	 *
	 * You know, what, I'm going to make this an iterator based so that I have
	 * a little bit of flexibility.
	 */
	template <class T_iter>
	void associateLoci(T_iter begin, const T_iter& end);

	/**
	 * Loading function - must be subclassed
	 * This function does the heavy lifting,region is not in the map and is the only function that needs
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


protected:
	// A map from id -> Region*
	unordered_map<uint,Region*> _region_map;
	// A map from alias -> id (to get region, essentially use
	// regionMap[aliasMap[alias]])
	unordered_map<string,set<Region*> > _alias_map;


	// OK, this is ugly!
	/**
	 * It's a map from chromosome -> map of intervals -> set of Regions
	 * Basically, use as the following:
	 * _region_bounds[chromosome][position] = {All Regions containing that position}
	 */
	unordered_map<short,interval_map<uint, set<Region*> > > _region_bounds;

private:
	/**
	 * Copy constructor + assignment operator.  This object should not be
	 * copied or assinged.  EVER!
	 */
	RegionCollection(const RegionCollection&);
	RegionCollection& operator=(const RegionCollection&);

	/**
	 * Deletes regions that have no SNPs associated with them
	 */
	void Squeeze();

	/**
	 * Special value used for if the requested region is not in the map
	 */
	Region region_not_found;
	/**
	 * Special value used for if the requested alias not present
	 * (NOTE: will be initialized by the default constructor)
	 */
	const set<Region*> empty_region_set;

};
}


#endif
