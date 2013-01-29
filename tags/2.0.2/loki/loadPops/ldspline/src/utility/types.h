//
// C++ Interface: types
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TYPES_H
#define TYPES_H

//Required for DRAND48
//#include <acconfig.h>

#include <stdint.h>
#include <vector>
#include <string>

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <stdlib.h>

#ifndef uint 
typedef unsigned int uint;
#endif

#ifdef WIN32

#define lrand48 lrand
#define srand48 srand
float drand48 ();
//float nanf(const char *tagp);



#include <float.h>
//Let's pretend MS isn't at least 8 years out of date!
//If you are stuck using MS' boken compiler, you will need these
//#define isnormal(x)  _fpclass(x)==_FPCLASS_NN || _fpclass(x)==_FPCLASS_PN
//#define isnan(x) _fpclass(x) == _FPCLASS_SNAN || _fpclass(x) == _FPCLASS_QNAN
#endif

#define MAX_LINE_LENGTH 600000
#ifdef TEST_APP
#include <gtest/gtest.h>
#endif


namespace LDUtility {



using namespace boost;

typedef std::vector<std::string> StringArray;

//typedef tokenizer<char_separator<char> > strtokenizer;

typedef dynamic_bitset<> BitSetType;
typedef std::vector<BitSetType> BitSetArray;

/**
 * @brief This is the storage medium for genetypes. 
 * It provides a label and bitarray
 */
struct GtStorage {
	std::string label;						///<This is used to identify a given genotype
	BitSetType individuals;				///<Which individuals in a given snp have this genotype
	GtStorage(const char *label, BitSetType individuals) : label(label), individuals(individuals) {}		
	GtStorage() : label("") {}
};
struct LocusID {
	uint chrID;
	uint locID;
	
	LocusID(uint chr, uint loc) : chrID(chr), locID(loc) { }
};

typedef std::vector<GtStorage> GenotypeArray;

//typedef uint uint32_t;

}

#endif

