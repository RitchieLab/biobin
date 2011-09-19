/* 
 * File:   genegenemodel.h
 * Author: torstees
 *
 * Stores gene indexes and implication index
 * We are storing the implication index with the indices, but assume that
 * there is a region manager present so that we can generate the actual models
 *
 * Created on March 8, 2011, 10:29 AM
 */

#ifndef GENEGENEMODEL_H
#define	GENEGENEMODEL_H

#include "snpsnpmodel.h"
#include <fstream>
#include "regionmanagerdb.h"
#include "snpdataset.h"
#include <soci.h>


namespace Knowledge {

class GeneGeneModel {
public:
	GeneGeneModel() : implicationIndex(0.0) { }
	GeneGeneModel(uint g1, uint g2, float implicationIndex);
	GeneGeneModel(const GeneGeneModel& orig) : implicationIndex(orig.implicationIndex), regions(orig.regions) { }
	virtual ~GeneGeneModel() {}

	/**
	 * For STL set, it actually sorts so largest II is first
    * @param other
    * @return true if II > other or II==other.II and region < other.region
    */
	bool operator<(const GeneGeneModel& other) const;
	bool operator==(const GeneGeneModel& other) const;
	GeneGeneModel& operator=(const GeneGeneModel& other);

	void Write(std::ostream& os, bool isBinary = false) const;
	void WriteBinary(std::ostream& os) const;

	bool Load(std::istream& file, bool isBinary = false);
	bool LoadBinary(std::istream& file, uint count = 2);

	std::string ToString(const char *sep = "\t")const;
	std::string ToString(const RegionManager& regions, const char *sep = "\t") const;

	uint GenerateModels(SnpSnpModel::Collection& models, RegionManager& regions) const;

	void GenerateRandomModels(uint count, SnpSnpModel::Collection& models, RegionManager& regions);

	uint EstimateModelCount(const RegionManager& regions) const;


	uint operator[](uint idx) const;				///< Access to the indexes
	uint Size();										///< Number of regions
	float implicationIndex;							///< Implication index

	friend std::ostream &operator<<(std::ostream &stream, GeneGeneModel model);
private:
	std::vector<uint> regions;						///< Indexes

};

inline
std::ostream& operator<<(std::ostream& out, GeneGeneModel model) {
	out<<model.ToString();
	return out;
}


inline
GeneGeneModel::GeneGeneModel(uint g1, uint g2, float ii) : implicationIndex(ii) {
	if (g1 == g2)
		regions.push_back(g1);
	else {
		if (g1 > g2) {
			uint t = g1;
			g1 = g2;
			g2 = t;
		}

		regions.push_back(g1);
		regions.push_back(g2);
	}
}

inline
GeneGeneModel& GeneGeneModel::operator=(const GeneGeneModel& other) {
	regions = other.regions;
	implicationIndex = other.implicationIndex;
	return *this;
}

inline
bool GeneGeneModel::operator<(const GeneGeneModel& other) const {

	//We actually want to sort these in the opposite direction (at least for II), since
	//we want the largest II on top of the set
	if (implicationIndex == other.implicationIndex)
		return regions < other.regions;
	else
		return implicationIndex > other.implicationIndex;
}

inline
bool GeneGeneModel::operator==(const GeneGeneModel& other) const {
	if (implicationIndex == other.implicationIndex)
		return regions == other.regions;
	return false;
}

inline
uint GeneGeneModel::operator[](uint idx) const {
	return regions[idx];
}

inline
uint GeneGeneModel::EstimateModelCount(const RegionManager& regions) const {
	return regions[this->regions[0]].SnpCount() * regions[this->regions[1]].SnpCount();
}

inline
uint GeneGeneModel::GenerateModels(SnpSnpModel::Collection& models, RegionManager& regions) const {
	Region& region = regions[this->regions[1]];
	return regions[this->regions[0]].GenerateModels(models, region, implicationIndex);
}
inline
void GeneGeneModel::GenerateRandomModels(uint count, SnpSnpModel::Collection& models, RegionManager& regions)  {
	Region& region = regions[this->regions[1]];
	regions[this->regions[0]].GenerateRandomModels(count, models, region, implicationIndex);
}

inline
uint GeneGeneModel::Size() {
	return regions.size();
}

inline
std::string GeneGeneModel::ToString(const char *sep) const {
	std::stringstream ss;
	ss<<Utility::Join(regions, "\t")<<"\t"<<implicationIndex;
	return ss.str();
}

inline
std::string GeneGeneModel::ToString(const RegionManager& regionMgr, const char *sep) const {
	std::stringstream ss;
	std::vector<uint>::const_iterator itr = regions.begin();
	std::vector<uint>::const_iterator end = regions.end();
	while (itr != end) {
		Region r = regionMgr[*itr++];
		r.WriteToArchive(ss, " ");
		ss<<"\t";
	}
	ss<<implicationIndex;
	return ss.str();
}

inline
void GeneGeneModel::Write(std::ostream& os, bool isBinary) const {
	if (isBinary)
		WriteBinary(os);
	else {
		os<<Utility::Join(regions, "\t")<<"\t"<<implicationIndex<<"\n";
	}
}


inline
bool GeneGeneModel::Load(std::istream& os, bool isBinary) {
	bool success = false;
	if (isBinary)
		return LoadBinary(os);
	else {
		static char line[65536];
		memset(line, '\0', 65536);
		os.getline(line, 65535);


		if (strlen(line) > 0) {
			Utility::StringArray words = Utility::Split(line, "\t");
			uint count = words.size();
			regions.clear();
			for (uint i=0; i<count-1; i++)
				regions.push_back(atoi(words[i].c_str()));
			implicationIndex = atof(words[count-1].c_str());
			success = true;
		}
	}
	return success;
}

inline
void GeneGeneModel::WriteBinary(std::ostream& os) const {
	uint v=0;
	for (uint i=0; i<regions.size(); i++) {
		v = regions[i];
		os.write((char*)&v, 4);
	}
	os.write((char*)&implicationIndex, 4);
}

inline
bool GeneGeneModel::LoadBinary(std::istream& os, uint count) {
	uint v=0;
	regions.clear();
	for (uint i=0; i<count; i++) {
		os.read((char*)&v, 4);
		regions.push_back(v);
	}
	os.read((char*)&implicationIndex, 4);
	return true;
}

}

#endif	/* GENEGENEMODEL_H */

