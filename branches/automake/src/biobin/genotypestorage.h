/* 
 * File:   genotypestorage.h
 * Author: torstees
 *
 * Created on June 21, 2011, 11:07 AM
 */

#ifndef GENOTYPESTORAGE_H
#define	GENOTYPESTORAGE_H

#include <map>
#include <sstream>
#include "utility/types.h"
#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace BioBin {
	

#define MISSING (uint)-1;	
	
/**
 * For normal SNPs, well store a 1 in each strand where there is a major allele.
 * For more complex genotypes (where there can be more than 2 alleles), we'll 
 * store the more complex genotype integer value in the map, complexGenotypes
 * 
 * The complex genotypes should all have an index > than bitsets.size()
 * 
 */
struct GenotypeStorage {
	GenotypeStorage();
	GenotypeStorage(uint genotypeCount);
	GenotypeStorage(const GenotypeStorage& other);
	
	void SetGenotype(uint index, char s1, char s2);
	void SetGenotype(uint index, int genotype);
	int GetGenotype(uint index) const;
	
	
	/// all bitsets should be sized == MaxGenotypeIndex
	boost::dynamic_bitset<> strand1;
	boost::dynamic_bitset<> strand2;
	boost::dynamic_bitset<> missing;
		
	/// index -> genotype
	std::map<uint, char> complexGenotypes;
	
	
	std::string GetGenotypes(const char *) const;
	
	/**
	 * This should be set up ahead of time...
	 */
	/// how many irregular alleles there are
	static std::vector<char> alleleCount;			

	uint GenotypeCount() const;
};


inline
GenotypeStorage::GenotypeStorage() 
		: strand1(0), 
		  strand2(0), 
		  missing(0) { }

inline
GenotypeStorage::GenotypeStorage(uint genotypeCount) 
		: strand1(genotypeCount), 
		  strand2(genotypeCount), 
		  missing(genotypeCount) {
}

inline
GenotypeStorage::GenotypeStorage(const GenotypeStorage& other)  {
	complexGenotypes = other.complexGenotypes;
	strand1 = other.strand1;
	strand2 = other.strand2;
	missing = other.missing;
}

inline
uint GenotypeStorage::GenotypeCount() const {
	return strand1.size();
}

inline
std::string GenotypeStorage::GetGenotypes(const char *sep) const {
	uint snpCount = GenotypeCount();
	std::stringstream ss;
	for (uint i=0; i<snpCount; i++) 
		ss<<GetGenotype(i)<<sep;
	
	return ss.str();
}
}
#endif	/* GENOTYPESTORAGE_H */

