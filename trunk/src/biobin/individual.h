/* 
 * File:   individual.h
 * Author: torstees
 *
 * Created on June 7, 2011, 12:53 PM
 */

#ifndef INDIVIDUAL_H
#define	INDIVIDUAL_H

#include <iostream>
#include <vector>
#include <string>
#include <map>

//#include "genotypestorage.h"

using std::string;
using std::map;

#include "knowledge/Locus.h"

namespace BioBin {

class Individual{

public:
	Individual(const string& name);
	
private:
	map<Knowledge::Locus*, short> _genotype_map;
};

/*
struct Individual {
public:
	Individual()																			
			: indID(""), pedID(""), status(-1) { }
	Individual(const std::string& indID, uint genotypeCount, uint binCount, const std::string& pedID = "")	
			: indID(indID), 
			pedID(pedID), 
			status(-1),
			genotypes(genotypeCount),
			binData(binCount, 0) {}
	Individual(const Individual& orig)												
			: indID(orig.indID), 
			pedID(orig.pedID), 
			status(-1),
			genotypes(orig.genotypes),
			binData(orig.binData) { }
	virtual ~Individual() {} 
	void Init(const std::string& indID, uint genotypeCount, uint binCount, const std::string& pedID = "") {
		this->indID = indID;
		genotypes = GenotypeStorage(genotypeCount);
		binData = std::vector<uint>(binCount, 0);
	}
	
	std::ostream& WriteGenotypes(std::ostream& file, const char *sep = " ") const;
	std::ostream& WriteBins(std::ostream& file, const char *sep=" ") const;
	uint BinCount(uint index) const;
	
	void ApplyBinCounts(std::vector<uint>& binCounts) {
		std::vector<uint>::iterator itr = binData.begin();
		std::vector<uint>::iterator end = binData.end();
		//std::cerr<<"ID: "<<indID<<" "<<Utility::Join(binData, " ")<<"\n";
		uint i=0;
		while (itr != end) 
			binCounts[i++]+= *itr++;
		//std::cerr<<"--"<<Utility::Join(binCounts, " ")<<"\n";
	}

	uint GenotypeCount() const { return genotypes.GenotypeCount(); }
	
	void IncrementBin(uint index, uint amount) {
		assert(index < binData.size());
		if ((int)amount > 0)
			binData[index]+=amount;
	}
	
	void SetGenotype(uint index, char value) {
		genotypes.SetGenotype(index, value);
	}
	

	 * Allow us to avoid sticking .0s at the end of MDR style status (which might
	 * mess up MDR...not sure)
	 /
	static bool ConvertStatusToInteger;			

	std::string indID;					///< Individual's ID from dataset
	std::string pedID;					///< Just in case...
	float status;							///< 
protected:
	GenotypeStorage genotypes;			///< This is where we keep all of the genotype data
	std::vector<uint> binData;			///< bin hits
};*/
/*
inline
std::ostream& Individual::WriteBins(std::ostream& file, const char *sep) const {
	if (ConvertStatusToInteger)
		file<<indID<<sep<<(int)status<<sep<<Utility::Join(binData, sep);
	else
		file<<indID<<sep<<status<<sep<<Utility::Join(binData, sep);
	return file;
}

inline
std::ostream& Individual::WriteGenotypes(std::ostream& file, const char *sep) const {
	if (ConvertStatusToInteger)
		file<<indID<<sep<<(int)status<<sep<<genotypes.GetGenotypes(sep);
	else
		file<<indID<<sep<<status<<sep<<genotypes.GetGenotypes(sep);
	return file;
}

inline
uint Individual::BinCount(uint index) const {
	return binData[index];
}
*/

}
#endif	/* INDIVIDUAL_H */

