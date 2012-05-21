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

/*!
 * \brief A class to manage a group of Regions.
 * This class manages a group of Regions by using the AddRegion method.  Users
 * should use this method fo creating Regions as opposed to the public
 * constructors of the Region class.  By using this class, you have a single
 * point of reference that a user can access a Region by ID or by alias
 * (including canonical name).  Ideally, there should be only one collection
 * of Regions, though there may come a time when multiple collections make sense.
 *
 * This class is an abstract base class so that implementation-specific loading
 * behavior can be defined.  For now, RegionCollectionSQLite is the only subclass
 * defined.
 *
 * This class is a complete rewrite of the RegionManager class in Biofilter 1.0
 * It is designed to eliminate the id->index confusion.
 */
class RegionCollection{
	
public:
	/*!
	 * Typedef to hide implementation of the collection of Regions.
	 */
	typedef set<Region*>::const_iterator const_region_iterator;

	/*!
	 * Create a new RegionCollection and initialize the special "region not found"
	 * object that is returned when needed.
	 */
	RegionCollection():region_not_found("Not Found",-1,-1, 0, 0){};
	/*!
	 * Destroy the RegionCollection object.  This will in turn also delete any
	 * Regions that were created with the AddRegion methods.
	 */
	virtual ~RegionCollection();

	/*!
	 * \brief Adds a Region to the collection.
	 * Adds a region to the collection where the true and effective start and
	 * stop values are identical.
	 *
     * \param name The primary name
     * \param id The id from the database (must be unique)
     * \param start	The true and effective start of the region to add
     * \param stop The true and effective stop of the region to add
	 * \param aliases Comma separated list of aliases
     */
	Region* AddRegion(const string& name, uint id, short chrom, uint start, uint stop, const string& aliases = "");

	/*!
	 * \brief Adds a region to the collection.
	 * Adds a Region to the collection where the true and effective start and
	 * stop values are different.
	 *
     * \param name The primary name
     * \param id The id from the database (must be unique)
     * \param effStart The effective start of the Region to add
     * \param effStop The effective stop of the Region to add
     * \param trueStart The true start of the Region to add
     * \param trueStop The true stop of the Region to add
	 * \param aliases Comma separated list of aliases
     */
	Region* AddRegion(const string& name, uint id, short chrom, uint effStart, uint effStop, uint trueStart, uint trueStop, const string& aliases = "");
	
	/*!
	 * \brief Access bracket operator, indexing by id.
	 * Returns a reference to a Region object that has the given ID.
	 * NOTE: contrary to the default behavior of a map, this will not insert
	 * values into the map.  To insert, you must use AddRegion.
	 *
	 * \param id The unique database ID of the region.
	 *
	 * \return A reference to the Region.
	 */
	Region& operator[](const uint id);

	/**
	 * Access bracket operator, indexing by alias
	 * NOTE: contrary to the default behavior of a map, this will not insert
	 * values into the map.  To insert, you must use AddRegion
	 */
	//const set<Region&>& operator[](const string&);

	/*!
	 * \brief Access bracket operator, indexing by id.
	 * Returns a reference to a Region object that has the given ID.  This
	 * version is for const-object access (typically lvalue access).
	 *
	 * \param id The unique database ID of the region.
	 *
	 * \return A reference to the Region.
	 */	
	const Region& operator[](const uint id) const;

	/*!
	 * \brief Access all Regions that have the given alias.
	 * Access by alias is controlled by an iterator, because sometimes aliases
	 * are not unique.  This method returns an iterator to the beginning of the
	 * set of Regions that are identified by the given alias.  If no Regions are
	 * identified by the given alias, the iterator is over the special empty
	 * set.
	 *
	 * \param alias The alias for the Region objects.
	 *
	 * \return An iterator to the start of the set of Regions identified by the
	 * given alias.
	 */
	const_region_iterator aliasBegin(const string& alias) const;
	/*!
	 * \brief Returns an iterator to the end of the set of Regions for an alias.
	 * This method returns the special one-past-the-end iterator for the Regions
	 * identified by the given alias.
	 *
	 * \param alias The alias for region objects.
	 *
	 * \return The end iterator for the collection of Regions identified by alias.
	 */
	const_region_iterator aliasEnd(const string& alias) const;

	/*!
	 * \brief Returns an iterator to a set of Regions at a given position.
	 * This function returns an iterator to the set of Regions that are located
	 * at a given position.  This can be used to find which Regions contain a
	 * given Locus, as Regions may overlap.  If no Regions are at the given
	 * location, the iterator is over the special empty set.
	 *
	 * \param chrom The chromosome (index) of the location
	 * \param pos The position in question.
	 *
	 * \return An iterator to the beginning of the collection of Regions
	 */
	const_region_iterator positionBegin(short chrom, uint pos) const;
	/*!
	 * \brief Returns an iterator to the end of a set of Regions at a position,
	 * Because Regions may overlap, this function returns the special end
	 * iterator of the set of Regions that contain a given position.
	 *
	 * \param chrom The chromosome (index) of the location
	 * \param pos The position in question.
	 *
	 * \return An iterator to the end of the collection of Regions
	 */
	const_region_iterator positionEnd(short chrom, uint pos) const;

	/*!
	 * \brief Determine if a given alias is valid.
	 * Determines if a given alias refers to any valid Regions.  Semantically
	 * equivalent to aliasBegin(alias) == aliasEnd(alias), but somewhat quicker.
	 *
	 * \param alias The alias to check for validity.
	 *
	 * \return A boolean that is true <==> alias is contained in the list of aliases.
	 */
	bool isValid(const string& alias) const {
		return _alias_map.find(alias) != _alias_map.end();
	}

	/*!
	 * \brief Compares the region agains the private "invalid" Region.
	 * This funciton is to be used to determine if a Region that was returned
	 * by an operator[] is a valid Region, or simply a proxy for a "not found"
	 * object.  This is needed b/c operator[] MUST return a Region object, so
	 * if the given ID did not match a Region in the collection, we would either
	 * have to throw an exception (costly and dangerous), or return a "Not Found"
	 * proxy object.  We have chosen the latter, and since the not found proxy
	 * object is private, we need to test for equality through a public function.
	 *
	 * \param other The Region to test for validity.
	 *
	 * \return A boolean that is true <==> other has the same ID as the special
	 * "not found" Region.
	 */
	bool isValid(const Region& other) const;

	/*!
	 * \brief Associates a collection of Loci to the Regions that contain them.
	 * Takes a list of Locus objects (realistically, just chrom->bp locations)
	 * and adds them to the Regions that contain them.  The template parameters
	 * must be incrementable iterators that return a Locus* object when
	 * dereferenced.
	 *
	 * \tparam begin A forward iterator of Locus* objects
	 * \tparam end The end iterator of a collection of Locus* objects
	 */
	template <class T_iter>
	void associateLoci(T_iter begin, const T_iter& end);

	/*!
	 * \brief The Loading function - must be subclassed.
	 * This function does the heavy lifting, loading a collection of regions
	 * that are identified by the given IDs or aliases.  If both are empty, load
	 * all Regions contained in the database.
	 *
	 * \param ids A list of IDs to load
	 * \param aliasList A list of aliases to include
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	virtual uint Load(const unordered_set<uint>& ids,
			const vector<string>& aliasList) = 0;

	/*!
	 * \brief Calls Load(...) with an empty ID list.
	 * Calls the Load function with an empty ID list.  This should not be
	 * subclassed unless specific behavior is desired.
	 *
	 * \param aliasList A list of aliases to include
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	uint Load(const vector<string>& aliasList);

	/*!
	 * \brief Calls Load(...) with an empty aliasList.
	 * Calls the Load function with an empty aliasList.  This should not be
	 * subclassed unless specific behavior is desired.
	 *
	 * \param ids A list of IDs to include
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	uint Load(const unordered_set<uint>& ids);

	/**
	 * Calls Load(...) with an empty set of ids and and empty aliasList
	 */
	//virtual uint Load(const uint popID);

	/*!
	 * \brief Calls Load(...) with empty set of ids.
	 * This function calls the Load function with an empty set of IDs (which
	 * should in turn call Load with an empty set of IDs and aliases).  Using
	 * this function should load ALL Regions from storage.
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	uint Load();


	//! The configured population to use here (default "NO-LD")
	static string pop_str;
	//! The amount of gene boundary expansion (if using no-ld, default 0)
	static int gene_expansion;

protected:
	//! A map from id -> Region*
	unordered_map<uint,Region*> _region_map;
	//! A map from alias -> set of Region*
	unordered_map<string,set<Region*> > _alias_map;


	// OK, this is ugly!
	/*!
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

	/*!
	 * Deletes regions that have no SNPs associated with them
	 */
	void Squeeze();

	/*!
	 * Special value used for if the requested region is not in the map
	 */
	Region region_not_found;

	/*!
	 * Special value used for if the requested alias not present
	 * (NOTE: will be initialized by the default constructor)
	 */
	const set<Region*> empty_region_set;

};

template <class T_iter>
void RegionCollection::associateLoci(T_iter begin, const T_iter& end){
	while (begin != end){
		interval_map<uint, set<Region*> >& chrom_map =
				_region_bounds[(*begin)->getChrom()];
		interval_map<uint, set<Region*> >::const_iterator region_itr =
				chrom_map.find((*begin)->getPos());
		if (region_itr != chrom_map.end()){
			set<Region*>::iterator set_itr=region_itr->second.begin();
			set<Region*>::const_iterator set_end=region_itr->second.end();
			while (set_itr != set_end){
				(*set_itr)->addLocus(*begin);
				++set_itr;
			}
		}
		++begin;
	}
}
}


#endif
