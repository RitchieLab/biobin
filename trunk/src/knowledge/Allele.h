/*
 * Allele.h
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_ALLELE_H
#define KNOWLEDGE_ALLELE_H

#include <string>
#include <set>
#include <ostream>
#include <stdlib.h>



namespace Knowledge{

/*!
 * \brief A single allele.
 * This class represents a single allele, containing the data, frequency, and
 * the position of the allele.  In this case, the position refers not to the
 * position on the genome, but rather which allele (1st, 2nd, etc.).
 * The position is typically 0-based with 0 being the major allele.  Clients of
 * this class must ensure that for a given locus, the position and data for
 * an allele are unique.
 */
class Allele{
public:

	/*!
	 * Constructs a new allele given the data, frequency and position.
	 * \param data The data associated with the allele (usually a single character, but not guaranteed)
	 * \param freq The frequency of this allele
	 * \param pos The position of the allele in the given input file
	 */
	Allele(const std::string& data, float freq, unsigned int pos) :
			_data(&(*(s_string_pool.insert(data).first))), _freq(freq), _pos(pos) {
	}

	/*!
	 * Destructor for the allele.  Currently a noop.
	 */
	~Allele() {}

	/*!
	 * Comparison operator, used in STL constructs.  Alleles are ordered by
	 * frequency, then by data.  This should ensure that for a given locus, the
	 * alleles are strictly ordered.  The comparison of alleles between loci is
	 * nonsensical.
	 *
	 * NOTE: The comparison of frequencies is reversed from standard, so the
	 * allele with the highest frequency will be first in an STL ordered
	 * container.  That is, other._freq > this._freq => other < this.
	 *
	 * \param other The allele to compare to
	 */
	bool operator<(const Allele& other) const;
	//bool operator>(const Allele&) const;
	//bool operator==(const Allele&) const;

	/*!
	 * Returns the frequency associated with this allele.
	 * \return The frequency of the allele.
	 */
	float getFreq() const { return _freq;}
	// Will we need to access the allele?  we'll see
	/*!
	 * Returns the data of the allele.
	 * \return The data associated with this specific allele.
	 */
	const std::string& getData() const { return *_data;}
	// Returns the position of the allele (0 is reference typically)
	unsigned short getPos() const {return _pos;}

	/*!
	 * \brief Prints the allele.
	 * This function prints the allele to the given ostream, using the
	 * separator given.  The format of the output is <data><sep><freq>.
	 * \param o The output stream to print to
	 * \param sep The separator used to separate data from frequency
	 */
	void print(std::ostream& o, const std::string& sep=":") const;

private:

	const std::string* _data;
	float _freq;
	unsigned short _pos;

	static std::set<std::string> s_string_pool;
};

}

/*!
 * \brief A shortcut for printing the allele.
 * This function allows for
 */
std::ostream& operator<<(std::ostream& o, const Knowledge::Allele& l);

#endif /* ALLELE_H_ */
