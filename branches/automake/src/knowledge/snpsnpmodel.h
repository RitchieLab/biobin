//
// C++ Implementation: snpsnpmodel.h
//
// Description: Defines functionality associated with Snp/Snp models
//
// Unlike much of the other knowledge data, these actually point to 
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SNP_SNP_MODEL_H
#define __SNP_SNP_MODEL_H
#include "biomodel.h"
#include <set>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <sstream>
#include "def.h"

namespace Knowledge {



class SnpSnpModel : public BioModel {
public:
	typedef std::set<SnpSnpModel> Collection;

	SnpSnpModel();
	SnpSnpModel(uint snp1, uint snp2, float implicationIndex);
	virtual ~SnpSnpModel();

	/* file buffer requirements */
	bool operator<(const SnpSnpModel& other) const;
	bool operator==(const SnpSnpModel& other) const;
	SnpSnpModel& operator=(const SnpSnpModel& other);
	void Load(std::istream& os, uint snpcount = 2);
	void Write(std::ostream& os) const;
	uint operator[](const uint idx) const;
	uint Size() const;
	std::string ToString(const char *sep = "\t") const;
	std::string ToFLString() const;
protected:
	void WriteBinary(std::ostream& file) const;
	bool LoadBinary(std::istream& file, uint modelSize = 2);
	std::vector<uint> snps;
};

inline
uint SnpSnpModel::Size() const {
	return snps.size();
}

inline
uint SnpSnpModel::operator[](const uint idx) const {
	return snps[idx];
}

inline
SnpSnpModel::SnpSnpModel() : BioModel(0.0) { }
inline
SnpSnpModel::SnpSnpModel(uint snp1, uint snp2, float implicationIndex) : BioModel(implicationIndex) {
	if (snp1 == snp2) {
		snps.push_back(snp1);
	}
	else	{
		if (snp1 > snp2) {
			uint t = snp1;
			snp1 = snp2;
			snp2 = t;
		}
		snps.push_back(snp1);
		snps.push_back(snp2);
	}
}

inline
SnpSnpModel::~SnpSnpModel() { }

inline
bool SnpSnpModel::operator<(const SnpSnpModel& other) const {
	return snps < other.snps;
}


inline
bool SnpSnpModel::operator==(const SnpSnpModel& other) const {
	return snps == other.snps;
}

inline
SnpSnpModel& SnpSnpModel::operator=(const SnpSnpModel& other) {
	implicationIndex = other.implicationIndex;
	snps = other.snps;
	return *this;
}

inline
std::ostream& operator<<(std::ostream& out, const SnpSnpModel& other) {
	other.Write(out);
	return out;
}

inline
std::string SnpSnpModel::ToString(const char *sep) const {
	std::stringstream ss;

	std::vector<uint>::const_iterator itr = snps.begin();
	std::vector<uint>::const_iterator end = snps.end();
	while (itr != end)
		ss<<*itr++<<sep;
	ss<<std::setprecision(1)<<ImplicationIndex();
	return ss.str();
}

inline
std::string SnpSnpModel::ToFLString() const {
	std::stringstream ss;

	std::vector<uint>::const_iterator itr = snps.begin();
	std::vector<uint>::const_iterator end = snps.end();
	while (itr != end)
		ss<<std::setw(MAX_RS_LENGTH)<<std::right<<*itr++;
	ss<<std::setw(MAX_II_LENGTH)<<std::right<<std::setprecision(1)<<ImplicationIndex();
	return ss.str();
}

inline
void SnpSnpModel::Write(std::ostream& os) const {
	if (BinaryArchive) {
		WriteBinary(os);
		return;
	}

	os<<ToFLString()<<"\n";
}

inline
void SnpSnpModel::WriteBinary(std::ostream& file) const {
	//uint modelSize = snps.size();
	std::vector<uint>::const_iterator itr = snps.begin();
	std::vector<uint>::const_iterator end = snps.end();
	
	while (itr != end) {
		uint locus = *itr++;
		file.write((char*)&locus, sizeof(uint));
	}
	file.write((char*)&implicationIndex, sizeof(float));           //# implication index
}

inline
void SnpSnpModel::Load(std::istream& os, uint modelSize) {
	snps.clear();

	if (BinaryArchive) {
		LoadBinary(os, modelSize);
		return;
	}



	for (uint i=0; i<modelSize; i++) {
		char rs[] = "            ";
		os.read(rs, MAX_RS_LENGTH);
		snps.push_back(atoi(rs));
	}
	char ii[] = "       ";
	os.read(ii, MAX_II_LENGTH);
	implicationIndex = atof(ii);

	//Get rid of the return character
	os.read(ii, 1);
}
inline
bool SnpSnpModel::LoadBinary(std::istream& file, uint modelSize) {
	//file.read((char*)&modelSize, sizeof(uint));

	for (size_t i=0; i<modelSize; i++) {
		uint locus = 0;
		file.read((char*)&locus, sizeof(uint));
		snps.push_back(locus);
	}
	file.read((char*)&implicationIndex, sizeof(float));
	return modelSize > 0;
}


}

#endif //__SNP_SNP_MODEL_H
