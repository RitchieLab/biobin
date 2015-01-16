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
#include <boost/pool/pool.hpp>
#include <new>

#include "Allele.h"
#include "LocusPostion.h"

namespace Knowledge {

/*!
 * \brief A class representing a locus in the genome.
 * This class represents a locus, which we have defined as a single position
 * on a chromosome that has some variation between individuals.  Traditionally,
 * these have been SNPs, but a Locus is somewhat more general, as we may have
 * any variant at this location.  This class will keep track of all the
 * different alleles possible at this position.
 *
 * NOTE: The user of this class must ensure that no two loci will have the
 * same chromosome and position.
 */
class Locus {

public:

	static void* operator new(size_t size);
	static void operator delete(void* deadObj, size_t size);

	//typedef std::vector<Allele>::const_iterator const_allele_iterator;

	/*!
	 * \brief Constructs a Locus object using the chromosome string.
	 * This constructor uses the chromosome string, along with the position to
	 * define a single locus.  Optionally, a user can provide an ID, which could
	 * correspond to an RSID.  If not provided, a unique ID will be created in
	 * the form "chr#-pos"
	 *
	 * \param chrom_str The chromosome string ID of the chromosome (see the _chrom_list variable for acceptable formatting)
	 * \param pos The position on the chromosome of this locus
	 * \param id The unique ID string of this locus (optional)
	 */
	Locus(const std::string& chrom_str, uint pos, const std::string& id = "");

	/*!
	 * \brief Constructs a locus object using the chromosome index.
	 * This constructor uses the chromosome index, which is what is directly
	 * stored in the object.  Typically, this is chromosome # - 1, where X=23,
	 * Y=24, XY=25, MT=26 (though this is subject to change, and will be
	 * different for species other than humans).  Again, an optional ID can be
	 * provided, and if not provided, one will be created.
	 *
	 * \param chrom The chromosome index of the chromosome
	 * \param pos The position on the chromosome of the locus
	 * \param id The unique ID string of this locus (optional)
	 */
	Locus(short chrom, uint pos, const std::string& id = "");

	/*!
	 * Return the ID of this Locus (passed in or auto-generated).
	 *
	 * \return The unique ID of the Locus.
	 */
	const std::string& getID() const {
		return _id;
	}
	;

	/*!
	 * Return the string identifying the chromosome (see _chrom_list)
	 *
	 * \return The chromosome string of the Locus.
	 */
	const std::string& getChromStr() const {
		return getChromStr(_chrpos.getChrom());
	}
	;

	/*!
	 * Return the chromosome index of this Locus (helpful for indexing)
	 *
	 * \return The cromosome index of the Locus.
	 */
	unsigned short getChrom() const {
		return _chrpos.getChrom();
	}
	;

	/*!
	 * Return the position of this Locus
	 *
	 * \return The base pair location of this Locus.
	 */
	unsigned int getPos() const { return _chrpos.getPos(); }

	/*!
	 * \brief Returns the distance to another Locus.
	 * Gives the distance (absolute value of the difference of the positions)
	 * from this Locus to another on the same chromosome.  If on different
	 * chromosomes, returns -1.
	 *
	 * \param other A Locus to determine distance
	 *
	 * \return The distance from this locus to the other (-1 if on different chromosomes)
	 */
	unsigned int distance(const Locus& other) const;

	/*!
	 * \brief Comparison operator for use in STL ordered classes.
	 * This compares 2 loci and orders them by chromosome then position.
	 *
	 * \param other The other Locus in the comparison
	 *
	 * \return
	 */
	bool operator<(const Locus& other) const;

	/*!
	 * \brief A function to print a Locus.
	 * This function prints a
	 */
	void print(std::ostream& o, const std::string& sep = ",") const;

	/*!
	 * \brief Converts a chromosome index into a chromosome string
	 * Helper function for converting a chromosome index to a chromosome string
	 * NOTE: This is the inverse of getChrom below.
	 *
	 * \param chrom A chromosome index
	 *
	 * \return A chromosome string
	 */
	static const std::string& getChromStr(unsigned short chrom);

	/*!
	 * \brief Converts chromosome string into an index
	 * Helper function for converting a chromosome string to an index.  Inverse
	 * of getChromStr
	 *
	 * \param chrom_str A chromosome string
	 * \return The corresponding index for the chromosome
	 */
	static unsigned short getChrom(const std::string& chrom_str);

	// Special flag for invalid chromosome (bad string or bad position)
	static const std::string invalid_chrom;

private:
	// No copying or assigning - use pointers, please!
	Locus(const Locus&);
	Locus& operator=(const Locus&);

	/*!
	 * \brief Helper method for creating an ID
	 * Method that creates an ID based on the chromosome and position and saves
	 * the created Id into the _id member.
	 */
	void createID();

	LocusPostion _chrpos;

	// index into list of chromosomes
	//short _chrom;
	// Position on the chromosome
	//unsigned int _pos;
	// Identifier of this Locus (could be a RSID or anything)
	std::string _id;

	// Vector of a list of chromosomes
	static const std::vector<std::string> _chrom_list;

	// memory pool of Locus objects (for speed, Locus objects are fairly lightweight)
	static boost::pool<> s_locus_pool;

};
/*
template<class Allele_itr>
void Locus::addAlleles(Allele_itr begin, const Allele_itr& end) {
	while (begin != end) {
		_alleles.push_back(*begin);
		++begin;
	}
}

template<class Str_itr>
void Locus::addAllelesStr(Str_itr begin, const Str_itr& end){
	unsigned short i=0;
	while (begin != end) {
		_alleles.push_back(Allele(*begin, i++));
		++begin;
	}
}
*/

}

std::ostream& operator<<(std::ostream& o, const Knowledge::Locus& l);

namespace std {

template<>
struct less<Knowledge::Locus*> {

	bool operator()(const Knowledge::Locus* x,
			const Knowledge::Locus* y) const {
		return (y != 0 && x != 0) ? (*x) < (*y) : x < y;
	}
};

template<>
struct less<const Knowledge::Locus*> {

	bool operator()(const Knowledge::Locus* x,
			const Knowledge::Locus* y) const {
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};

}

#endif /* KNOWLEDGE_LOCUS_H */
