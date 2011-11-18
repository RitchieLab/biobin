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

using std::string;
using std::set;
using boost::unordered_map;

class Locus;
class Group;

namespace Knowledge{

class Region{

public:

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
	 * Equality operator
	 */
	bool operator==(const Region& other) const;

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
	unordered_map<uint, set<const Group*> > _group_map;

	string _name;
	short _chrom;
	uint _id;
	uint _true_start;
	uint _true_end;
	uint _eff_start;
	uint _eff_end;

};
}


#endif /* KNOWLEDGE_REGION_H */
