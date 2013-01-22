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

	typedef std::vector<Allele>::const_iterator const_allele_iterator;

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
	Locus(const std::string& chrom_str, uint pos, bool rare = false,
			const std::string& id = "");

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
	Locus(short chrom, uint pos, bool rare = false, const std::string& id = "");

	/*!
	 * \brief Add an associated allele to this particular Locus.
	 * This function adds an allele to this locus, given the data and the
	 * frequency.  The position of the allele (see Allele) is autoincremented,
	 * starting at 0.  It is incumbent upon the user to ensure that the
	 * frequency of all alleles will sum to 1.
	 *
	 * NOTE: It is assumed that the user will add all alleles before accessing
	 * any of the encoding or decoding of genotype information.
	 *
	 * \param allele The data of the given allele
	 * \param freq The frequency of the allele to add.
	 */
	void addAllele(const std::string& allele, float freq);

	/*!
	 * \brief Adds a list of alleles to this locus.
	 * This function adds a list of alleles to this locus.  This is used under
	 * the assumption that this list is coming from another Locus
	 *
	 * \param begin An iterator pointing to the beginning of the allele list
	 * \param end An iterator pointing to the end of the allele list
	 */
	template<class Allele_itr>
	void addAlleles(Allele_itr begin, const Allele_itr& end);

	/*!
	 * Returns an iterator to the beginning of the list of alleles.  Used when
	 * copying Loci, especially in a liftOver situation.
	 * \return An iterator pointing to the beginning of the allele List
	 */
	const_allele_iterator beginAlleles() const {
		return _alleles.begin();
	}
	/*!
	 * Returns an iterator to the end of the list of alleles.  Used when
	 * copying Loci, especially in a liftOver situation.
	 * \return An iterator pointing to the end of the allele List
	 */
	const_allele_iterator endAlleles() const {
		return _alleles.end();
	}

	/*!
	 * Sets the major allele, which may be different from the calculated major
	 * allele (in the case of overall-major-allele)
	 * \param majAllele The string of the major allele
	 */
	void setMajorAllele(const std::string& majAllele);

	/*!
	 * \brief Return the major allele frequency (that which is most common).
	 * Returns the most common allele among all those that have been entered
	 * already.  This is determined by the ordering of the Allele object.
	 *
	 * \return the major allele frequency.
	 */
	float majorAlleleFreq() const;

	/*!
	 * \brief Return the minor allele frequency (sum of all other alleles)
	 * Returns the minor allele frequency.  In the case of a biallelic Locus,
	 * this is simply the frequency of the allele that is present in the
	 * minority of the population.  In the case of the multi-allelic genes, this
	 * is the sum of all other alleles other than the major allele, or
	 * 1-majorAlleleFreq().  Note that it is possible that the minor allele
	 * frequency could be greater than 1/2 in the case of multi-allelic genes.
	 *
	 * \return The minor allele frequency.
	 */
	float minorAlleleFreq() const;

	/*!
	 * \brief Returns the rarity status of the Locus.
	 * This returns whether or not this locus is considered rare accoring to the
	 * input parameter during creation of the Locus.  Note that even Loci with
	 * a high minor allele frequency can be considered rare because in the case
	 * population, the minor allele frequency could be low.
	 */
	bool isRare() const {
		return _is_rare;
	}

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
		return getChromStr(_chrom);
	}
	;

	/*!
	 * Return the chromosome index of this Locus (helpful for indexing)
	 *
	 * \return The cromosome index of the Locus.
	 */
	short getChrom() const {
		return _chrom;
	}
	;

	/*!
	 * Return the position of this Locus
	 *
	 * \return The base pair location of this Locus.
	 */
	unsigned int getPos() const {
		return _pos;
	}
	;

	/*!
	 * \brief Returns the position of the major allele.
	 * Returns the position (alternate #) of the major allele.  Usually, this
	 * will be 0 (reference allele), but this is not a guarantee.  Typically,
	 * will be used in determining if a bin contains a minor allele "hit".
	 *
	 * \return The alternate number of the major allele.
	 */
	unsigned short getMajorPos() const {
		return (*_alleles.begin()).getPos();
	}

	/*!
	 * \brief Determine if the given allele is a minor allele or not
	 * This function checks the given allele data against the major allele to
	 * determine if the given allele is the major allele or not.  Note that if
	 * the given allele is not actually in the allele list, this function will
	 * still return "true", as it is not the same as the major allele.
	 *
	 * \param allele The allele data to check against the major allele
	 *
	 * \return A boolean that is false <==> the given allele is the major allele
	 */
	bool isMinor(const std::string& allele) const;

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
	 * \brief Encodes the genotype in a single value
	 * This function encodes a given genotype, given as the positions a1 and a2,
	 * which can the be decoded at a later time.  This is encoded as (# of
	 * alleles)*a1 + a2 (similar to a bitmask).
	 *
	 * NOTE: After encoding a genotype, no further alleles should be added, or
	 * the results will be incorrect.
	 *
	 * \param a1 The allele # on the 1st strand
	 * \param a2 The allele # on the 2nd strand
	 *
	 * \return The encoded genotype
	 */
	short encodeGenotype(unsigned int a1, unsigned int a2) const;

	/*!
	 * \brief Decodes the genotype into a pair of values.
	 * This function takes an encoded genotype and decodes it into two separate
	 * allele positions.  This function is the inverse of encodeGenotype.
	 *
	 * NOTE: If alleles have been added after the genotype was encoded, this
	 * function will give incorrect results.
	 *
	 * \param encoded_type The encoded genotype from encodeGenotype
	 *
	 * \return A pair of integers representing the allele numbers.
	 */
	std::pair<unsigned int, unsigned int> decodeGenotype(short encoded_type) const;

	/*!
	 * \brief A function to print a Locus.
	 * This function prints a
	 */
	void print(std::ostream& o, const std::string& sep = ",",
			bool printAlleles = false) const;

	/*!
	 * \brief a function to print the alleles
	 */
	void printAlleles(std::ostream& o, const std::string& sep = "|") const;

	/*!
	 * \brief Converts a chromosome index into a chromosome string
	 * Helper function for converting a chromosome index to a chromosome string
	 * NOTE: This is the inverse of getChrom below.
	 *
	 * \param chrom A chromosome index
	 *
	 * \return A chromosome string
	 */
	static const std::string& getChromStr(short chrom);

	/*!
	 * \brief Converts chromosome string into an index
	 * Helper function for converting a chromosome string to an index.  Inverse
	 * of getChromStr
	 *
	 * \param chrom_str A chromosome string
	 * \return The corresponding index for the chromosome
	 */
	static short getChrom(const std::string& chrom_str);

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

	// index into list of chromosomes
	short _chrom;
	// Position on the chromosome
	unsigned int _pos;
	// Identifier of this Locus (could be a RSID or anything)
	std::string _id;

	// A set of alleles for this Locus
	// IN this case, the final element will be the largest, so
	// *(_alleles.rbegin()) is the major allele
	std::vector<Allele> _alleles;

	//flag determining rarity
	bool _is_rare;

	// Vector of a list of chromosomes
	static const std::vector<std::string> _chrom_list;

	// memory pool of Locus objects (for speed, Locus objects are fairly lightweight)
	static boost::pool<> s_locus_pool;

};

template<class Allele_itr>
void Locus::addAlleles(Allele_itr begin, const Allele_itr& end) {
	while (begin != end) {
		_alleles.push_back(*begin);
		++begin;
	}
}

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
