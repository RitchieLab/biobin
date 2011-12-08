/* 
 * File:   converter.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:25 PM
 */

#ifndef KNOWLEDGE_LIFTOVER_CONVERTER_H
#define	KNOWLEDGE_LIFTOVER_CONVERTER_H

#include <vector>
#include <string>
#include <ostream>
#include <map>

//#include "Chain.h"

using std::multimap;
using std::string;
using std::vector;

namespace Knowledge {

class Locus;

namespace Liftover{

class Chain;

using Knowledge::Locus;

class Converter {
public:
	Converter(const string& origBuild, const string& newBuild);
	virtual ~Converter();

	void addChain(const string& chainDetails);
	
	void ConvertDataset(const string& mapFilename, multimap<Locus*, Locus*>& converted);
	/**
	 * Takes a list (iterable) of Locus object ptrs, applies the chaining and
	 * returns the converted Locus objects in the supplied container.  Also
	 * logs any dropped Locus objects to the supplied ostream object.
	 */
	template <class T_iter, class T_cont>
	void ConvertDataset(T_iter itr, const T_iter& end, T_cont& newBuild, std::ostream& droppedLog);

	/**
	 * Takes a map filename and returns the new Locus objects in teh supplied
	 * container.
	 */
	template <class T_cont>
	void ConvertDataset(const string& mapFilename, T_cont& newBuild, std::ostream& droppedLog);

	/**
	 * Takes a list (iterable) of Locus object ptrs and chains them together and
	 * returns the result as a multimap
	 */
	template <class T_iter>
	void ConvertDataset(T_iter itr, const T_iter& end,
			multimap<Locus*, Locus*>& converted);

	/**
	 * Loads all of the chain files from a given source (pure virtual)
	 */
	virtual int Load() = 0;


protected:
	multimap<short, Chain*> _chains;					///< Chrom -> chains
	string _origBuild;
	string _newBuild;

private:
	// No copying or assignment allowed
	Converter(const Converter& orig);
	Converter& operator=(const Converter& other);

	// Function to read a map filename and return the results in a vector of
	// Locus Objects
	void readMapFile(const string& mapFilename, vector<Locus*>& array_out) const;
};


}
}

#endif	/* CONVERTER_H */

