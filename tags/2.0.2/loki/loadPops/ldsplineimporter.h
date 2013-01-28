/* 
 * File:   ldSplineImporter.h
 * Author: torstees
 *
 * Created on September 17, 2010, 1:50 PM
 */

#ifndef LDSPLINEIMPORTER_H
#define	LDSPLINEIMPORTER_H

#include "ldspline/ldspline.h"
using namespace Spline;

#include <sqlite3.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>

using std::string;
using std::map;
using std::vector;
using std::cerr;

class LdSplineImporter {
private:

	enum CutoffType { D_PRIME, R_SQUARED};

	struct PopulationSpline {
		string name;					///< CEU/JPT/etc
		string desc;					///< comment to help inform users who might not be familiar with the 3 letter names
		string filename;			///< The filename associated with the splines

		string GetPopulationName(CutoffType statType, float value) const {
			std::stringstream ss;
			ss<<name<<"-"<<(statType == R_SQUARED ? "RS" : (statType == D_PRIME ? "DP" : "UNK"))<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setprecision(2)<<value;
			return ss.str();
		}

		string GetPopulationName(const string& statType, float value) const{
			std::stringstream ss;
			ss<<name<<"-"<<statType<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setprecision(2)<<value;
			return ss.str();
		}
		PopulationSpline(std::string name, std::string desc, std::string filename) : name(name), desc(desc), filename(filename) {}
	};

	struct RegionBoundary {
		int geneID;
		int lower;
		int upper;
		int source_id;
		RegionBoundary(int geneID, int lower, int upper, int src) :
			geneID(geneID), lower(lower), upper(upper), source_id(src){}
	};


public:

	LdSplineImporter(const string& config_fn, const string& db_fn);

	LdSplineImporter(const string& config_fn, sqlite3 *db_conn);

	// Load the population given a DB connection
	void loadPops();

	~LdSplineImporter();


private:
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

	void UpdateZones();
	void ProcessLD(LocusLookup& chr, const PopulationSpline& sp, const map<string, int>& popIDs);
	void LoadGenes();
	void InitPopulationIDs(std::map<std::string, int>& popIDs, const PopulationSpline& sp);

	// Retrieves and returns the indexes on a table.  Also, drops the indexes
	void getAndDropIndexes(const string& tbl_name, map<string, string>& indexes_out);
	void restoreIndexes(const string& tbl_name, const map<string, string>& index_map);

	std::vector<PopulationSpline> splines;			///<population -> ldspline filename
	std::vector<std::pair<CutoffType, float> > cutoffs;
	//std::vector<float> dp;								///<The various DPrime values we are splining on
	//std::vector<float> rs;								///<The various RSquared values we are splining on

	// A mapping of chromosome -> vector or regions (so we only have to load them once!)
	std::map<short, vector<RegionBoundary> > _region_map;
	std::vector<RegionBoundary> regions;

	// Database connection
	sqlite3* _db;
	bool _self_open;
	bool _write_db;
	string dbFilename;

	//temp table name
	static const string _tmp_bnd_tbl;

	// Vector of a list of chromosomes
	static const vector<string> _chrom_list;

	// convert chromosome strings to DB equivalent
	static short getChrom(const string& chrom_str);

	// DB processing functions
	static int parseGenes(void*, int, char**, char**);
	static int parseSingleInt(void*, int, char**, char**);
	static int parseRegionIndex(void*, int, char**, char**);

};


#endif	/* LDSPLINEIMPORTER_H */

