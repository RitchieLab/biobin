/* 
 * File:   ldSplineImporter.h
 * Author: torstees
 *
 * Created on September 17, 2010, 1:50 PM
 */

#ifndef LDSPLINEIMPORTER_H
#define	LDSPLINEIMPORTER_H
#include <soci.h>
#include "ldspline/ldspline.h"
#include <sstream>
#include <iomanip>

using namespace Spline;

namespace Biofilter {

class LdSplineImporter {
public:

	struct PopulationSpline {
		std::string name;					///< CEU/JPT/etc
		std::string desc;					///< comment to help inform users who might not be familiar with the 3 letter names
		std::string filename;			///< The filename associated with the splines

		std::string GetPopulationName(const char *statType, float value){
			std::stringstream ss;
			ss<<name<<"-"<<statType<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setprecision(2)<<value;
			return ss.str();
		}
		PopulationSpline(std::string name, std::string desc, std::string filename) : name(name), desc(desc), filename(filename) {}
	};

	struct RegionBoundary {
		int geneID;
		int lower;
		int upper;
		std::string chrom;
		RegionBoundary(int geneID, std::string chrom, int lower, int upper) : geneID(geneID), lower(lower), upper(upper), chrom(chrom) {}
	};

	LdSplineImporter();					///< Construction

	virtual ~LdSplineImporter();
	/**
	 * @brief Parse configuration
    * @param filename
	 *
	 * Example:
	 * rs 0.9 0.8 0.6
	 * dp 0.9 0.8 0.6
	 * CEU /path/to/ceu.ldspline Descriptive note about CEU population
	 * JPT /path/to/jpg.ldspline Descriptive note about the population
	 * ...
    */
	void LoadConfiguration(const char *filename);
	void LoadGenes(soci::session& sociDB, const char *chrom);
	int InitPopulation(soci::session& sociDB, const char *pop, const char *popDesc);
	void InitPopulationIDs(soci::session& sociDB, std::map<std::string, int>& popIDs, PopulationSpline& sp, const char *type, std::vector<float>& stats);

	void ProcessLD(soci::session& sociDB, LocusLookup& chr, PopulationSpline& sp, std::map<std::string, int>& popIDs);
	void Process(soci::session& sociDB);
	int GetPopID(soci::session& sociDB, const char *type, float threshold);
private:
	std::vector<PopulationSpline> splines;			///<population -> ldspline filename
	std::vector<float> dp;								///<The various DPrime values we are splining on
	std::vector<float> rs;								///<The various RSquared values we are splining on
	std::vector<RegionBoundary> regions;

};


}

#endif	/* LDSPLINEIMPORTER_H */

