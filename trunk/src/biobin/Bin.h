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
	//! Typedef to hide implementation details of the Locus containment
	typedef set<Knowledge::Locus*>::const_iterator const_locus_iterator;

	/*!
	 * \brief Construct a bin representing a Group (or Pathway)
	 *
	 * \param pop_mgr The population manager containing information about
	 * the variants and people.
	 * \param grp The group that this Bin represents
	 */
	Bin(const PopulationManager& pop_mgr, Knowledge::Group* grp);
	/*!
	 * \brief Construct a Bin representing a Region (or Gene)
	 *
	 * \param pop_mgr The population manager containing information about
	 * the variants and people.
	 * \param reg The Region that this Bin represents
	 */
	Bin(const PopulationManager& pop_mgr, Knowledge::Region* reg);
	/*!
	 * \brief Construct an Intergenic Bin
	 *
	 * \param pop_mgr The population manager containing information about
	 * the variants and people.
	 * \param chrom The chromosome of the bin
	 * \param bin The bin number, counting from 0 from the beginning of the
	 * chromosome.
	 */
	Bin(const PopulationManager& pop_mgr, short chrom, int bin);

	/*!
	 * \brief Return the size of the bin.
	 *
	 * \return The size (number of contributions) of this bin.
	 */
	int getSize() const;
	/*!
	 * \brief Return the number of variants in the bin.
	 *
	 * \return The number of variants (Loci) that the bin contains
	 */
	int getVariantSize() const {return _variants.size();}
	/*!
	 * \brief Return a unique ID (for type) of the bin.
	 *
	 * \return The ID of the
	 */
	int getID() const{
		return _is_intergenic ? _member.bin_no :
				(_is_group ? _member.group->getID() : _member.region->getID());
	}

	Knowledge::Group* getGroup() const { return _is_group ? _member.group : NULL;}

	const string& getName() const{ return _name;}
	bool isGroup() const {return _is_group;}
	bool isIntergenic() const {return _is_intergenic;}
	short getChrom() const {return _chrom;}

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
