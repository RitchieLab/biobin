/*
 * Region.h
 *
 *  Created on: Nov 18, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_REGION_H
#define KNOWLEDGE_REGION_H

#include <string>
#include <set>
#include <map>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/iterator/iterator_facade.hpp>

using std::string;
using std::set;
using std::map;
using std::vector;
using boost::unordered_map;



namespace Knowledge{

class Locus;
class Group;

/*!
 * \brief A class representing a region on a chromosome.
 * This class represents a region on a chromosome, defined to be a segment on
 * a chromosome that has a defined beginning and ending position.  A region is
 * basically a generalized gene.  The beginning and end positions may have a
 * "true" and "effective" beginning and end, which can account for population
 * differences or gene extension.  The "true" boundaries are those defined in
 * the database for the typical population, while the "effective" boundaries
 * are those that are population-specific or extended.
 */
class Region{

public:

	// TODO: templatize this if needed
	/*!
	 * \brief An iterator over the groups
	 * This iterator will traverse either all groups of a given category, or
	 * all groups associated with the given Region.  Constructed via the
	 * groupBegin and groupEnd methods of the Region class.
	 *
	 * NOTE: This iterator only implements a forward traversal, so incrementing
	 * is the only oepration available.
	 */
	class const_group_iterator : public boost::iterator_facade<const_group_iterator, Group* const, boost::forward_traversal_tag>{
	public:
		//const_group_iterator() {}

		/*!
		 * \brief Create a const_group_iterator that iterates over all groups
		 * This is the constructor for iterating over all groups in the given
		 * group_map.  If begin is false, this is constructing an iterator that
		 * points to the end of the iteration.
		 *
		 * \param group_map A mapping of IDs to sets of groups to iterate over
		 * \param begin A flag telling if we want a start or end iterator
		 */
		const_group_iterator(const map<uint, set<Group*> >& group_map, bool begin) :
				_is_global(true){
			_map_iter = group_map.begin();
			_map_end = group_map.end();
			
			if(_map_end == _map_iter){
				//OK, we got passed an empty map!
				if(begin){
					_set_iter = _empty_set.begin();
				} else {
					_set_iter = _empty_set.end();
				}
				_curr_end = _empty_set.end();
			}else{
				// Not empty
				if (begin){
					_set_iter = (*_map_iter).second.begin();
					_curr_end = (*_map_iter).second.end();
				} else {
					_map_iter = _map_end;
					--_map_iter;
					_set_iter = (*_map_iter).second.end();
					_curr_end = (*_map_iter).second.end();
					++_map_iter;
				}
				
			}
		}

		/*!
		 * \brief Create a const_group_iterator for iterating over groups from a single source.
		 * This constructor creates a const_group_iterator for iterating over
		 * groups from a single source.  If the group_id is not in the group_map,
		 * we assume that the iterator iterates over an empty set.
		 *
		 * \param group_map The mapping of source IDs to sets of groups
		 * \param group_id The ID of the source of interest
		 * \param begin A flag determining if the iterator is a start or end iterator
		 */
		const_group_iterator(const map<uint, set<Group*> >& group_map, uint group_id, bool begin) : _is_global(false) {
			// We don't need the _map_iter here
			map<uint, set<Group*> >::const_iterator g_itr = group_map.find(group_id);
			if(g_itr == group_map.end()){
				if(begin){
					_set_iter = _empty_set.begin();
				}else{
					_set_iter = _empty_set.end();
				}
				_curr_end = _empty_set.end();
			}else{
				if (begin){
					_set_iter = (*g_itr).second.begin();
				}else{
					_set_iter = (*g_itr).second.end();
				}
				_curr_end = (*g_itr).second.end();
			}
		}

	private:
		friend class boost::iterator_core_access;

		/*!
		 * \brief The incrementing function for iterators
		 * This increments an iterator, moving on to the next group, if needed
		 */
		void increment() {
			++_set_iter;
			if(_is_global && _set_iter == _curr_end){
				++_map_iter;
				if(_map_iter != _map_end){
					_set_iter = (*_map_iter).second.begin();
					_curr_end = (*_map_iter).second.end();
				}
			}
		}

		bool equal(const const_group_iterator& other) const{
			return this->_set_iter == other._set_iter;
		}

		Group* const& dereference() const {return *_set_iter;}

		// This should ALWAYS be empty!
		static const set<Group*> _empty_set;

		set<Group*>::const_iterator _set_iter;
		set<Group*>::const_iterator _curr_end;
		map<uint, set<Group*> >::const_iterator _map_iter;
		map<uint, set<Group*> >::const_iterator _map_end;
		bool _is_global;


	};

	// If I need a non-const iterator, I'll use it here
	//typedef set<string>::iterator alias_iterator;
	/*!
	 * A typedef to abstract the implementation of the storage of aliases
	 */
	typedef set<string>::const_iterator const_alias_iterator;

private:

	class Boundary{

	private:
		Boundary(short chr, uint start, uint stop):
			_chr(chr), _start(start), _stop(stop) {}

	public:
		friend class Region;

		bool operator<(const Boundary& other) const {
			return _chr == other._chr ?
					(_start == other._start ?
							_stop < other._stop : _start < other._start) :
					_chr < other._chr;
		}

		short getChrom() const { return _chr; }

	private:
		short _chr;
		uint _start;
		uint _stop;
	};

public:

	Region(const string& name, uint id);

	/*!
	 * \brief Construct a region using ony start/stop values.
	 * Constructs a region in which the effective boundaries are the same as
	 * the true boundaries of the region.  Here, the chromosome must be indexed,
	 * and the id refers to the databse ID of the region, which must be unique
	 * among all regions.  To find the index of the chromosome, you can use the
	 * Locus class.  Additionally, the name must be unique among the entire
	 * application.
	 *
	 * \param name The canonical name of the region (aliases can be added elsewhere)
	 * \param id The unique ID of this region
	 * \param chrom The indexed chromosome on which the region lies
	 * \param start The true and effective beginning of the region
	 * \param end The true and effective ending of the region
	 */
	Region(const string& name, uint id, short chrom, uint start, uint end);

	/*!
	 * \brief Construct a region using start/stop and effective (population) start/stop values.
	 * Constructs a region with potentially different true and effective
	 * boundaries.  Again, the id must be unique and the chromosome is indexed
	 * as from the static Locus methods.
	 *
	 * \param name The unique canonical name of the region
	 * \param id The unique ID of the region
	 * \param chrom The indexed chromosome on which the region lies
	 * \param eff_start The effective start of the region
	 * \param eff_end The effective end of the region
	 * \param true_start The true start of the region
	 * \param true_end The true end of the region
	 */
	Region(const string& name, uint id, short chrom,
			uint eff_start, uint eff_end, uint true_start, uint true_end);

	/*!
	 * \brief Associate a single locus with this region.
	 * This function assocates the given locus with this region.  Loci are
	 * stored as a set of pointers to Locus objects.  Typically loci are
	 * associated with a region because they are contained within it, though
	 * this need not be true; a locus can be associated with a region because it
	 * somehow affects the region or is in high LD with the region.
	 *
	 * \param locus the Locus object to associate with this region
	 */
	void addLocus(const Locus* locus);

	void addPopulationBoundary(short chr, uint start, uint stop){
		_pop_bounds.push_back(Boundary(chr, start, stop));
	}

	void addDefaultBoundary(short chr, uint start, uint stop){
		_def_bounds.push_back(Boundary(chr, start, stop));
	}

	/*!
	 * \brief A mass association of many Locus objects.
	 * Adds a whole lot of Loci (using iterators) with this region.  The
	 * iterators must be forward iterators, where the end iterator is the C++
	 * standard "one-past-the-end" iterator.  When dereferenced, the iterators
	 * must return a Locus* object.
	 *
	 * \tparam start The beginning iterator of the collection to add
	 * \tparam end The end iterator of the collection to add
	 */
	template <class T_iter>
	void addLoci(T_iter start, const T_iter& end);

	/*!
	 * \brief Associate a group with this region.
	 * Associate a single group with this region.  Note, we must supply a group
	 * type in addition to the actual group itself; a group type is typically
	 * a means of tracking which source this relationship is derived from.
	 *
	 * \param type The type (or source) of the Group object
	 * \param container The group containing this region
	 */
	void addGroup(uint type, Group& container);

	/*!
	 * \brief A mass association of Groups of a single type with this region.
	 * Adds a whole lot of Groups (using iterators) with this region.  The
	 * iterators again must be forward traversal iterators, which when
	 * dereferenced return Group* objects.  Note that this is strictly for
	 * a single type or source of groups.
	 *
	 * \param type The type (or source) of the Groups
	 * \tparam start The beginning iterator of the container of groups
	 * \tparam end The ending iterator of the container of groups
	 */
	template <class T_iter>
	void addGroups(uint type, T_iter& start, const T_iter& end);

	/*!
	 * \brief returns a string of aliases associated with this region.
	 * Returns a string of aliases associated with this region, separated by
	 * the given separator, which is a comma by default.
	 *
	 * \param sep The separator used to separate different aliases
	 * \return A sep-delinieated string of all aliases (including the canonical
	 * name) associated with this region.
	 */
	string getAliasString(const string& sep=",") const;

	/*!
	 * \brief Adds a list of aliases to this region.
	 * Adds a list of aliases to this region, where the aliases are separated
	 * with the given separator string.  Typically, this string will come from
	 * a database query, and if sep is given as an empty string, it will add the
	 * input string as a single alias.
	 *
	 * \param aliases A sep-separated string of aliases for this region
	 * \param sep The string separator
	 */
	void addAliases(const string& aliases, const string& sep=",");

	void addAlias(const string& alias) { _aliases.insert(alias); }

	/*!
	 * Return the number of Loci associated with this region
	 *
	 * \return The size of the Locus map
	 */
	uint locusCount() const {return _locus_map.size();}

	/*!
	 * Get the chromosome of this region
	 *
	 * \return The chromosome index corresponding to this region
	 */
	short getChrom() const {return _def_bounds.size() == 0 ? -1 : _def_bounds[0].getChrom();}

	/*!
	 * Get the ID of this region (used by the database and for indexing)
	 *
	 * \return The unique ID of this region
	 */
	uint getID() const {return _id;}

	/*!
	 * Return the canonical name of the region
	 *
	 * \return The canonical name of the region.
	 */
	const string& getName() const {return _name;}

	/*!
	 * \brief Check if this region contains the given Locus.
	 * This function checks the given Locus to see if it is associated with this
	 * Region.  To be considered a match, the ID string of the Locus must be
	 * the same as one that was added with addLocus or addLoci.  Note that if
	 * the ID strings are truly unique, this causes no problems, but if the
	 * ID strings are not unique, this could lead to spurious results.  Also, a
	 * user must be consistent with the IDs.  Ideally, the Locus object passed
	 * in here would be one of the same Locus objects actually added (as Locus
	 * objects are non-copyable and non-assignable).
	 *
	 * \param other The Locus to check for containment.
	 *
	 * \return A boolean that is true if the given Locus is associated with
	 * this region.
	 */
	bool containsLocus(const Locus& other) const;

	/*!
	 * Return an iterator to the beginning of the alias list
	 *
	 * \return A const iterator pointing to the beginning of the alias list
	 */
	const_alias_iterator aliasBegin() const {return _aliases.begin();}
	/*!
	 * Return an iterator to the end of the alias list.
	 *
	 * \return A const iterator pointing to the end of the alias list.
	 */
	const_alias_iterator aliasEnd() const {return _aliases.end();}

	/*!
	 * \brief Return an iterator over all groups.
	 * Returns an iterator to the beginning that when incremented, will iterate
	 * over all groups associated with this region.  Must be used with groupEnd
	 * to check for completion
	 *
	 * \return An iterator over all groups
	 */
	const_group_iterator groupBegin() const {return const_group_iterator(_group_map,  true);}

	/*!
	 * \brief Return an iterator over groups from a particular source.
	 * Returns an iterator to the beginning of a set of groups from a single
	 * source.  Must be used in concert with the corresponding groupEnd
	 *
	 * \return An iterator over groups from a given source
	 */
	const_group_iterator groupBegin(uint group) const {return const_group_iterator(_group_map, group, true);}

	/*!
	 * \brief Returns the global end iterator.
	 * Returns the iterator that is one past the end of the global iterator.
	 *
	 * \return The iterator that signals the global end of iterating over all groups
	 */
	const_group_iterator groupEnd() const {return const_group_iterator(_group_map, false);}

	/*!
	 * \brief Returns the local end iterator.
	 * Returns an iterator that signals the end of a list of groups from a given
	 * source.
	 *
	 * \return The iterator that is the end of a set of groups from a source.
	 */
	const_group_iterator groupEnd(uint group) const {return const_group_iterator(_group_map, group, false);}

	/*!
	 * \brief Comparison operator, for any STL as needed.
	 * Compares two regions.  Regions are ordered by chromosome, then by
	 * (true) starting position, then by (true) ending position.  Thus, two
	 * regions with different effective boundaries but identical true boundaries
	 * are considered to be equal.
	 *
	 * \param other The Region to compare this to
	 *
	 * \return A boolean determining if this is less than the other
	 */
	bool operator<(const Region& other) const;


private:

	// No copying or assignment allowed!
	Region(const Region& other);
	Region& operator=(const Region& other);

	// A map of IDs -> Loci
	unordered_map<string, const Locus*> _locus_map;

	// A list of all aliases to this Region
	set<string> _aliases;

	// A list / mapping of all groups associated w/ this region
	map<uint, set<Group*> > _group_map;

	string _name;
	uint _id;

	vector<Boundary> _pop_bounds;
	vector<Boundary> _def_bounds;

};


template <class T_iter>
void Region::addLoci(T_iter begin, const T_iter& end){
	while (begin != end){
		addLocus(*(*begin));
		++begin;
	}
}

template <class T_iter>
void Region::addGroups(uint type, T_iter& begin, const T_iter& end){
	set<Group*>& group_set = _group_map[type];
	while(begin != end){
		group_set.insert(*begin);
		++begin;
	}
}


}

// Overloading < operator for region pointers
namespace std{

template<>
struct less<Knowledge::Region*> {

	bool operator()(const Knowledge::Region* x, const Knowledge::Region* y) const{
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};


}
#endif /* KNOWLEDGE_REGION_H */
