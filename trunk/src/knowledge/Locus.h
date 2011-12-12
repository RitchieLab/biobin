/*
 * Locus.h
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_LOCUS_H
#define KNOWLEDGE_LOCUS_H

#include <string>
#include <set>
#include <vector>
#include <utility>
#include <ostream>
#include <stdlib.h>

#include "Allele.h"

using std::string;
using std::set;
using std::vector;
using std::pair;
using std::ostream;

namespace Knowledge{

class Locus{

public:

	/**
	 * Constructs a Locus object using the chromosome string (preferrable)
	 */
	Locus(const string& chrom_str, uint pos, const string& id="");

	/**
	 * Constructs a locus object using the chromosome index.
	 * (This is typically chromosome # - 1, where X=23, Y=24, XY=25, MT=26
	 */
	Locus(short chrom, uint pos, const string& id="");

	/**
	 * Add an associated allele to this particular Locus
	 */
	void addAllele(const string& allele, float freq);

	/**
	 * Return the major allele frequency (that which is most common)
	 */
	float majorAlleleFreq() const;
	/**
	 * Return the minor allele frequency (sum of all other alleles)
	 */
	float minorAlleleFreq() const;

	/**
	 * Return the ID of this Locus (passed in)
	 */
	const string& getID() const {return _id;};
	/**
	 * Return the string identifying the chromosome (see _chrom_list)
	 */
	const string& getChromStr() const {return getChromStr(_chrom);};
	/**
	 * Return the chromosome index of this Locus (helpful for indexing)
	 */
	short getChrom() const {return _chrom;};
	/**
	 * Return the position of this Locus
	 */
	uint getPos() const {return _pos;};
	/**
	 * Returns the position (alternate #) of the major allele.  Usually, this
	 * will be 0 (reference allele), but not always
	 */
	uint getMajorPos() const {return (*_alleles.begin()).getPos();}

	/**
	 * Determine if the given allele is a minor allele or not
	 */
	bool isMinor(const string& allele) const;

	/**
	 * Gives the distance from this Locus to another on the same chromosome.
	 * If on different chromosomes, returns -1
	 */
	uint distance(const Locus& other) const;
	/**
	 * Comparison operator for use in STL ordered classes.
	 * Ordered by chromosome then position.
	 */
	bool operator<(const Locus& other) const;

	/**
	 * Encodes the genotype in a single value, which can the be decoded
	 */
	short encodeGenotype(uint a1, uint a2) const;

	/**
	 * Decodes the genotype in a pair of values.  Inverse of encodeGenotype
	 */
	pair<uint, uint> decodeGenotype(short encoded_type) const;

	void print(ostream& o, const string& sep=",", bool printAlleles=false) const;

	void printAlleles(ostream& o, const string& sep="|") const;

	/**
	 * Helper function for converting a chromosome index to a chromosome string
	 * NOTE: This is the inverse of getChrom below.
	 */
	static const string& getChromStr(short chrom);
	/**
	 * Helper function for converting a chromosome string to an index
	 * Inverse of getChromStr
	 */
	static short getChrom(const string& chrom_str);

	// Special flag for invalid chromosome (bad string or bad position)
	static const string invalid_chrom;

private:
	// No copying or assigning - use pointers, please!
	Locus(const Locus&);
	Locus& operator=(const Locus&);

	/**
	 * Method that creates an ID based on the chromosome and position and saves
	 * the created Id into the _id member.
	 */
	void createID();

	// index into list of chromosomes
	short _chrom;
	// Position on the chromosome
	uint _pos;
	// Identifier of this Locus (could be a RSID or anything)
	string _id;

	// A set of alleles for this Locus
	// NOTE: I'm using set here to maintain sorted order
	// IN this case, the final element will be the largest, so
	// *(_alleles.rbegin()) is the major allele
	set<Allele> _alleles;

	// Vector of a list of chromosomes
	static const vector<string> _chrom_list;
};

}

ostream& operator<<(ostream& o, const Knowledge::Locus& l);

namespace std{

template<>
struct less<Knowledge::Locus*> {

	bool operator()(const Knowledge::Locus* x, const Knowledge::Locus* y) const{
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};

template<>
struct less<const Knowledge::Locus*> {
	
	bool operator()(const Knowledge::Locus* x, const Knowledge::Locus* y) const{
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};

}

#endif /* KNOWLEDGE_LOCUS_H */
