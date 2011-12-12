/*
 * Bin.h
 *
 *  Created on: Dec 5, 2011
 *      Author: jrw32
 */

#ifndef BIOBIN_BIN_H
#define BIOBIN_BIN_H

#include <set>

#include "knowledge/Group.h"
#include "knowledge/Region.h"
#include "knowledge/Locus.h"

using std::set;

namespace BioBin{

class PopulationManager;

class Bin{

public:
	typedef set<Knowledge::Locus*>::const_iterator const_locus_iterator;

	Bin(const PopulationManager& pop_mgr, Knowledge::Group*);
	Bin(const PopulationManager& pop_mgr, Knowledge::Region*);
	Bin(const PopulationManager& pop_mgr, short chrom, int bin);

	int getSize() const;
	int getVariantSize() const {return _variants.size();}
	int getID() const{
		return _is_intergenic ? -1 :
				(_is_group ? _member.group->getID() : _member.region->getID());
	}
	const string& getName() const{ return _name;}

	void addLocus(Knowledge::Locus* to_ins) {_variants.insert(to_ins);}

	/**
	 * Strict ordering is given by Groups, then Regions, then Intergenic, which
	 * are then sorted (by ID, chrom/position, chrom/position)
	 */
	bool operator<(const Bin& other) const;

	const_locus_iterator variantBegin() const {return _variants.begin();}
	const_locus_iterator variantEnd() const {return _variants.end();}

private:
	union{
		Knowledge::Group* group;
		Knowledge::Region* region;
		int bin_no;
	} _member;

	bool _is_group;
	bool _is_intergenic;
	short _chrom;

	mutable int _size_cache;
	mutable bool _cached;

	string _name;

	set<Knowledge::Locus*> _variants;

	const PopulationManager& _pop_mgr;
};
}

// Overloaded < for pointers
namespace std{

template<>
struct less<BioBin::Bin*>{
	bool operator() (const BioBin::Bin* x, const BioBin::Bin* y) const{
		return (*x) < (*y);
	}
};

}

#endif /* BIN_H_ */
