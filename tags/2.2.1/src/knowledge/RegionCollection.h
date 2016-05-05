#ifndef KNOWLEDGE_REGIONCOLLECTION_H
#define KNOWLEDGE_REGIONCOLLECTION_H

#include <string>
#include <vector>
#include <set>

#include "any_iterator.hpp"

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "Region.h"

namespace Knowledge{

class Information;

class Locus;

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
	
protected:

	//NOTE: The following 2 classes are designed to allow for the creation
	//of a RegionCollection from an arbitrary collection of forward-traversable
	//Locus* objects
	class Container {

	public:
		typedef IteratorTypeErasure::any_iterator<Knowledge::Locus* const, boost::forward_traversal_tag> const_iterator;

		virtual ~Container() {}
		virtual const_iterator begin() const = 0;
		virtual const_iterator end() const = 0;
	};

	template <class T_cont>
	class LocusContainer : public Container {
	public:

		LocusContainer(const T_cont& c) : _data(c) {}
		virtual ~LocusContainer() {}

		virtual const_iterator begin() const {return static_cast<const_iterator>(_data.begin());}
		virtual const_iterator end() const {return static_cast<const_iterator>(_data.end());}

	private:
		const T_cont& _data;
	};

public:
	/*!
	 * A means to iterate over the elements in the collection consistent with
	 * being in a set
	 */
	class const_iterator : public boost::iterator_facade<const_iterator, Region* const, boost::forward_traversal_tag>{

	public:
		const_iterator(boost::unordered_map<unsigned int, Region*>::const_iterator itr) : _itr(itr){}

	private:
		friend class boost::iterator_core_access;

		void increment() { ++_itr;}
		bool equal(const const_iterator& other) const { return _itr == other._itr;}
		Region* const& dereference() const { return (*_itr).second;}

		boost::unordered_map<unsigned int, Region*>::const_iterator _itr;

	};

	/*!
	 * Typedef to hide implementation of the collection of Regions.
	 */
	typedef std::set<Region*>::const_iterator const_region_iterator;

	/*!
	 * Create a new RegionCollection and initialize the special "region not found"
	 * object that is returned when needed.
	 */
	template<class T_cont>
	RegionCollection(const T_cont& loci) :
			_dataset(new LocusContainer<T_cont>(loci)),
			region_not_found("Not Found", -1) {
	}

	/*!
	 * Destroy the RegionCollection object.  This will in turn also delete any
	 * Regions that were created with the AddRegion methods.
	 */
	virtual ~RegionCollection();


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
	Region& operator[](const unsigned int id);

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
	const Region& operator[](const unsigned int id) const;

	/*!
	 * \brief An iterator to the beginning of the collection of Regions.
	 * This is an iterator the the beginning of the collection of Regions, but
	 * the iterator hides the fact that the implementation of the collection
	 * is a map (ie, *begin() == (*_region_map.begin()).second)
	 *
	 * \return an iterator to the beginning of the collection.
	 */
	const_iterator begin() const {return const_iterator(_region_map.begin());}

	/*!
	 * \brief An iterator to the end of the collection of Regions.
	 * This is an iterator the the end of the collection of Regions, but
	 * the iterator hides the fact that the implementation of the collection
	 * is a map (ie, *end() == (*_region_map.end()).second)
	 *
	 * \return an iterator to the end of the collection.
	 */
	const_iterator end() const {return const_iterator(_region_map.end());}

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
	const_region_iterator aliasBegin(const std::string& alias) const;
	/*!
	 * \brief Returns an iterator to the end of the set of Regions for an alias.
	 * This method returns the special one-past-the-end iterator for the Regions
	 * identified by the given alias.
	 *
	 * \param alias The alias for region objects.
	 *
	 * \return The end iterator for the collection of Regions identified by alias.
	 */
	const_region_iterator aliasEnd(const std::string& alias) const;

	const_region_iterator locusBegin(const Locus* loc) const;
	const_region_iterator locusEnd(const Locus* loc) const;

	/*!
	 * \brief Determine if a given alias is valid.
	 * Determines if a given alias refers to any valid Regions.  Semantically
	 * equivalent to aliasBegin(alias) == aliasEnd(alias), but somewhat quicker.
	 *
	 * \param alias The alias to check for validity.
	 *
	 * \return A boolean that is true <==> alias is contained in the list of aliases.
	 */
	bool isValid(const std::string& alias) const {
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
	virtual unsigned int Load(const boost::unordered_set<unsigned int>& ids,
			const std::vector<std::string>& aliasList) = 0;

	virtual void loadFiles() = 0;

	/*!
	 * \brief Calls Load(...) with an empty ID list.
	 * Calls the Load function with an empty ID list.  This should not be
	 * subclassed unless specific behavior is desired.
	 *
	 * \param aliasList A list of aliases to include
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	unsigned int Load(const std::vector<std::string>& aliasList);

	/*!
	 * \brief Calls Load(...) with an empty aliasList.
	 * Calls the Load function with an empty aliasList.  This should not be
	 * subclassed unless specific behavior is desired.
	 *
	 * \param ids A list of IDs to include
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	unsigned int Load(const boost::unordered_set<unsigned int>& ids);

	/*!
	 * \brief Calls Load(...) with empty set of ids.
	 * This function calls the Load function with an empty set of IDs (which
	 * should in turn call Load with an empty set of IDs and aliases).  Using
	 * this function should load ALL Regions from storage.
	 *
	 * \return 0 if successfully loaded, anything else upon error.
	 */
	unsigned int Load();


	//! The configured population to use here (default "NO-LD")
	static std::string pop_str;
	//! The amount of gene boundary expansion (if using no-ld, default 0)
	static int gene_expansion;
	//! A vector of strngs containing custom regions
	static std::vector<std::string> c_region_files;

	//! A vector of regions to filter on
	static std::vector<std::string> c_region_filter;

protected:
	/*!
	 * Inserts the given region into the id and alias maps
	 */
	void insertRegion(Region& region);

	//void insertRegionBound(Region& region, short chr, uint start, uint end);

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
	 *
	 * \return The region added (or the Region found by the id)
     */
	Region* AddRegion(const std::string& name, unsigned int id, short chrom,
			unsigned int start, unsigned int stop, const std::string& aliases = "");

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
	 *
	 * \return The region added (or the region found by the id)
     */
	Region* AddRegion(const std::string& name, unsigned int id, short chrom,
			unsigned int effStart, unsigned int effStop, unsigned int trueStart,
			unsigned int trueStop, const std::string& aliases = "");

	//! A map from id -> Region*
	boost::unordered_map<unsigned int,Region*> _region_map;
	//! A map from alias -> set of Region*
	boost::unordered_map<std::string,std::set<Region*> > _alias_map;

	//! object to get generalized information
	Information* _info;

	const Container* const _dataset;

	// Instead of mapping by position, let's map by locus!
	boost::unordered_map<const Locus*, std::set<Region*> > _locus_map;

private:
	/**
	 * Copy constructor + assignment operator.  This object should not be
	 * copied or assinged.  EVER!
	 */
	RegionCollection(const RegionCollection&);
	RegionCollection& operator=(const RegionCollection&);

	/*!
	 * Special value used for if the requested region is not in the map
	 */
	Region region_not_found;

	/*!
	 * Special value used for if the requested alias not present
	 * (NOTE: will be initialized by the default constructor)
	 */
	const std::set<Region*> empty_region_set;

};

}


#endif
