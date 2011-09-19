/* 
 * File:   genotypestorage.cpp
 * Author: torstees
 * 
 * Created on June 21, 2011, 11:07 AM
 */

#include "genotypestorage.h"

namespace BioBin {
std::vector<char> GenotypeStorage::alleleCount;
// s1,s2 are 0 if they have the major allele
void GenotypeStorage::SetGenotype(uint index, char s1, char s2) {
	if (s1 == (char)-1 || s2 == (char)-1 )
		missing[index] = true;
	else {
		strand1[index] = s1 > 0;
		strand2[index] = s2 > 0;

		//We'll only record complex genotypes that use anything other than the two
		//most common alleles
		if (s1 > 1 || s2 > 1) {
			complexGenotypes[index] = s2 + s1*alleleCount[index];
		}

	}
}

/**
 * By using the binary representation of the two bitstrings,
 * we can store only those genotypes with 1 or both major alleles,
 * and the case where we have two of the second most frequent allele.
 * This means that we will only add to complex genotypes whenever someone
 * has at least 1 of the really rare alleles-which will hopefully be only
 * a few per individual
 * @param index
 * @param genotype
 */
void GenotypeStorage::SetGenotype(uint index, int genotype) {
	if (index >= strand1.size()) {
		std::cerr<<" Oooops! Index is too big. WTF? "<<index<<"\t"<<strand1.size()<<"\n";
	}
	assert(index < strand1.size());
	if (genotype == -1)
		missing[index] = true;
	else {
		if (genotype == 0)
			strand1[index] = strand2[index] = false;
		else if (genotype == 1) {
			strand1[index] = 0;
			strand2[index] = 1;
		} else if (genotype == alleleCount[index]) {
			strand1[index] = 1;
			strand2[index] = 0;
		} else if (genotype == alleleCount[index] + 1) {
			strand1[index] = 1;
			strand2[index] = 1;
		} else
			complexGenotypes[index] = genotype;
	}
}

int GenotypeStorage::GetGenotype(uint index) const {
	if (missing[index])
		return MISSING;
	std::map<uint, char>::const_iterator item = complexGenotypes.find(index);
	std::map<uint, char>::const_iterator end = complexGenotypes.end();

	if (item != end)
		return item->second;
	return strand1[index]*alleleCount[index] + strand2[index];
}

}


