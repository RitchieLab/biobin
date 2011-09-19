/* 
 * File:   genegenemodelarchive.h
 * Author: torstees
 *
 * Created on March 11, 2011, 1:29 PM
 */

#ifndef GENEGENEMODELARCHIVE_H
#define	GENEGENEMODELARCHIVE_H

#include "genegenemodel.h"
#include "regionmanager.h"
#include "def.h"
#include <map>
#include <vector>
#include <fstream>

namespace Knowledge {

class GeneGeneModelArchive {
public:
	typedef std::set<GeneGeneModel>::iterator iterator;
	GeneGeneModelArchive() {}
	GeneGeneModelArchive(const GeneGeneModelArchive& orig);
	virtual ~GeneGeneModelArchive() {}

	/**
	 * Constructs the models associated with a given set of regionIDs
    * @param ids
    * @param regions
    */
	void AddRegions(Utility::IdCollection& ids, RegionManager& regions);


	/**
	 * Builds basis for model report
	 */
	uint SummarizeModelCounts(std::map<float, uint>& scores, const RegionManager& regions);
	uint Size();
	iterator Begin();
	iterator End();

	void Reset();

	void LoadFromArchive(std::istream& os, bool useBinary);
	void LoadFromArchive(const char *filename, bool useBinary);
	void WriteToArchive(std::ostream& os, bool useBinary);
	void WriteToArchive(const char *filename, bool useBinary);
	void WriteToArchive(const char *filename, const RegionManager& regions, std::map<float, uint>& scores, bool useBinary);
	void GenerateModels(SnpSnpModel::Collection& snpBasedModels, RegionManager& regions);
	static uint minImplicationIndex;				///< Used for writing only the "best" models to archives
	static uint maxModelCount;						///< Also used to restrict the size of the SNP/SNP archive
private:
	
	std::set<GeneGeneModel> models;
};


inline
uint GeneGeneModelArchive::Size() {
	return models.size();
}

inline
GeneGeneModelArchive::GeneGeneModelArchive(const GeneGeneModelArchive& other) {
	models = other.models;
}

inline
GeneGeneModelArchive::iterator GeneGeneModelArchive::Begin() {
	return models.begin();
}

inline
GeneGeneModelArchive::iterator GeneGeneModelArchive::End() {
	return models.end();
}

inline
void GeneGeneModelArchive::Reset() {
	models.clear();
}



}

#endif	/* GENEGENEMODELARCHIVE_H */

