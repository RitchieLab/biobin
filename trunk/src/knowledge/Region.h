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
#include <boost/unordered_map.hpp>
#include <boost/iterator/iterator_facade.hpp>

using std::string;
using std::set;
using boost::unordered_map;



namespace Knowledge{

class Locus;
class Group;

class Region{

public:

	// TODO: templatize this if needed
	class const_group_iterator : public boost::iterator_facade<const_group_iterator, Group* const, boost::forward_traversal_tag>{
	public:
		const_group_iterator() {}
		const_group_iterator(const unordered_map<uint, set<Group*> >::const_iterator& group_itr,
				const unordered_map<uint, set<Group*> >::const_iterator& map_end) : _is_global(true){
			_map_iter = group_itr;
			_map_end = map_end;
			if(_map_iter == _map_end){
				// In this case, we don't iterate over anything!
				_set_iter = _empty_set.begin();
				_curr_end = _empty_set.end();
			}else{
				_set_iter = (*_map_iter).second.begin();
				_curr_end = (*_map_iter).second.end();
			}
		}
		const_group_iterator(const unordered_map<uint, set<Group*> >& group_map, uint group_id, bool begin) : _is_global(false) {
			// We don't need the _map_iter here
			unordered_map<uint, set<Group*> >::const_iterator g_itr = group_map.find(group_id);
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

		bool equal(const const_group_iterator& other){
			return this->_set_iter == other._set_iter;
		}

		Group* const dereference() const {return *_set_iter;}

		// This should ALWAYS be empty!
		static const set<Group*> _empty_set;

		set<Group*>::const_iterator _set_iter;
		set<Group*>::const_iterator _curr_end;
		unordered_map<uint, set<Group*> >::const_iterator _map_iter;
		unordered_map<uint, set<Group*> >::const_iterator _map_end;
		unordered_map<uint, set<Group*> >::const_iterator _look_ahead;
		bool _is_global;


	};

	// If I need a non-const iterator, I'll use it here
	//typedef set<string>::iterator alias_iterator;

	typedef set<string>::const_iterator const_alias_iterator;

	/**
	 * Construct a region using ony start/stop values.
	 */
	Region(const string& name, uint id, short chrom, uint start, uint end);
	/**
	 * Construct a region using start/stop and effective (population) start/stop values
	 */
	Region(const string& name, uint id, short chrom,
			uint eff_start, uint eff_end, uint true_start, uint true_end);

	/**
	 * Associate a single locus with this region
	 */
	void addLocus(const Locus& locus);

	/**
	 * Add a whole lot of Loci (using iterators) with this region
	 * NOTe: We assume the iterators are iterating over a "set" of Locus* objects
	 */
	template <class T_iter>
	void addLoci(T_iter& start, const T_iter& end);

	/**
	 * Associate a single group with this region.  Note, we must supply a group
	 * type in addition to the actual group itself.
	 */
	void addGroup(uint type, Group& container);

	/**
	 * Add a whole lot of Groups (using iterators) with this region
	 * NOTe: We assume the iterators are iterating over a "set" of Group* objects
	 */
	template <class T_iter>
	void addGroups(uint type, T_iter& start, const T_iter& end);

	/**
	 * Return a string of aliases associated with this region, separated by
	 * the given separator
	 */
	string getAliasString(const string& sep=",") const;

	/**
	 * Add a list of aliases to this region, where the aliases are separated
	 * with the given separator string
	 */
	void addAliases(const string& aliases, const string& sep=",");

	/**
	 * Return the number of Loci associated with this region
	 */
	uint locusCount() const {return _locus_map.size();}

	/**
	 * Get the chromosome of this region
	 */
	short getChrom() const {return _chrom;}
	/**
	 * Get the ID of this region (used by the database and for indexing)
	 */
	uint getID() const {return _id;}
	/**
	 * Get the true start of the region
	 */
	uint getTrueStart() const {return _true_start;}
	/**
	 * Get the true stop of the region
	 */
	uint getTrueEnd() const {return _true_end;}
	/**const
	 * Get the effective (population) start of the region
	 */
	uint getEffStart() const {return _eff_start;}
	/**
	 * Get the effective (population) ending of the region
	 */
	uint getEffEnd() const {return _eff_end;}

	const string& getName() const {return _name;}

	/**
	 * Check if this region contains the given Locus (IDs must match)
	 */
	bool containsLocus(const Locus& other) const;

	/**
	 * Iterators over aliases
	 */
	const_alias_iterator aliasBegin() const {return _aliases.begin();}
	const_alias_iterator aliasEnd() const {return _aliases.end();}

	/**
	 * Iterators over groups
	 */
	const_group_iterator groupBegin() const {return const_group_iterator(_group_map.begin(), _group_map.end());}
	const_group_iterator groupBegin(uint group) const {return const_group_iterator(_group_map, group, true);}
	const_group_iterator groupEnd() const {return const_group_iterator(_group_map.end(), _group_map.end());}
	const_group_iterator groupEnd(uint group) const {return const_group_iterator(_group_map, group, false);}

	/**
	 * Comparison operator, for any STL as needed
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
	unordered_map<uint, set<Group*> > _group_map;

	string _name;
	short _chrom;
	uint _id;
	uint _true_start;
	uint _true_end;
	uint _eff_start;
	uint _eff_end;

};
}

// Overloading < operator for region pointers
namespace std{

template<>
struct less<Knowledge::Region*> {

	bool operator()(const Knowledge::Region* x, const Knowledge::Region* y){
		return (*x) < (*y);
	}
};

}
#endif /* KNOWLEDGE_REGION_H */
