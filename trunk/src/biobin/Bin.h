/*
 * Bin.h
 *
 *  Created on: Dec 5, 2011
 *      Author: jrw32
 */

#ifndef BIOBIN_BIN_H
#define BIOBIN_BIN_H

#include <set>
#include <list>

#include "knowledge/Group.h"
#include "knowledge/Region.h"
#include "knowledge/Locus.h"



namespace BioBin{

class PopulationManager;

class Bin{

public:
	//! Typedef to hide implementation details of the Locus containment
	typedef std::list<Knowledge::Locus*>::const_iterator const_locus_iterator;
	typedef std::list<Knowledge::Locus*>::iterator locus_iterator;

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
	 * \brief Construct a copy of a given bin.
	 * This constructs a copy of a bin, but IT DOES NOT COPY THE VARIANTS!!
	 * This is intended for use when constructing sub-bins based on a current
	 * bin (i.e., breaking out into sub-gene or functional information).  This
	 * is intentionally explicit in order to prevent implicitly generated copies
	 * from being created.
	 *
	 * Essentially, this is a convenience function, and you sould never be
	 * copying Bins.
	 *
	 * \param other The other bin to base this on.
	 */
	explicit Bin(const Bin& other);

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
	Knowledge::Region* getRegion() const {return (!_is_group && !_is_intergenic) ? _member.region : NULL; }

	const string& getName() const{ return _name;}
	bool isGroup() const {return _is_group;}
	bool isIntergenic() const {return _is_intergenic;}
	short getChrom() const {return _chrom;}

	void addLocus(Knowledge::Locus* to_ins) {_variants.push_back(to_ins);}

	/**
	 * Strict ordering is given by Groups, then Regions, then Intergenic, which
	 * are then sorted (by ID, chrom/position, chrom/position)
	 *
	 * Note: if the group / region / intergenic position are identical, then
	 * the operator looks at the _name member.
	 */
	bool operator<(const Bin& other) const;

	const_locus_iterator variantBegin() const {return _variants.begin();}
	const_locus_iterator variantEnd() const {return _variants.end();}
	locus_iterator variantBegin() { return _variants.begin();}
	locus_iterator variantEnd() {return _variants.end();}

	/*!
	 * \brief removes a Locus* from the bin.
	 * This function exactly mirrors the erase functionality of the STL set,
	 * which means that the itr will be non-functional after the erase.
	 *
	 * \param itr The iterator pointing to the element to erase
	 */
	locus_iterator erase(locus_iterator itr) {_cached = false; return _variants.erase(itr); }

	/*!
	 * \brief Adds extraconversion data to the bin.
	 * This function adds extra data to the bin.  The "extra data" is really
	 * just a string that we tack on to the end of the identifier.  Note that
	 * we'll need to change the name here, too
	 */
	void addExtraData(const string& addl_data) { _extra_data.push_back(addl_data); _name += addl_data;}


private:
	// No assignment, please!
	Bin& operator=(const Bin& other);

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

	std::string _name;
	std::vector<std::string> _extra_data;

	std::list<Knowledge::Locus*> _variants;

	const PopulationManager& _pop_mgr;
};
}

// Overloaded < for pointers
namespace std{

template<>
struct less<BioBin::Bin*>{
	bool operator() (const BioBin::Bin* x, const BioBin::Bin* y) const{
		return (y != 0 && x != 0) ? (*x) < (*y) : x < y;
	}
};

}

#endif /* BIN_H_ */
