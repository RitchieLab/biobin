/* 
 * File:   ldsplineimporter.cpp
 * Author: torstees
 * 
 * Created on September 17, 2010, 1:50 PM
 */
#include "ldsplineimporter.h"
#include <fstream>
#include "ldspline/ldspline.h"
using namespace std;
using namespace soci;


namespace Biofilter {


LdSplineImporter::LdSplineImporter() {
}

LdSplineImporter::~LdSplineImporter() {
}
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
void LdSplineImporter::LoadConfiguration(const char *filename) {
	std::ifstream file(filename);
	while (file.good() && !file.eof()) {
		char line[4096];
		file.getline(line, 4096);

		std::stringstream ss(line);

		std::istream_iterator<std::string> itr(ss);

		std::vector<std::string> tokens(itr, std::istream_iterator<std::string>());

		if (tokens.size() > 0) {
			if (tokens[0] == "rs" || tokens[0] == "RS") {
				std::vector<std::string>::iterator values = tokens.begin();
				std::vector<std::string>::iterator tokenEnd = tokens.end();
				while (++values != tokenEnd) {
					rs.push_back(atof(values->c_str()));
					cerr<<"RS: "<<*values<<"\t"<<rs[rs.size()-1]<<"\n";

				}
				cerr<< rs.size() << "Total RS values to be used.";
			}
			else if (tokens[0] == "dp" || tokens[0] == "DP") {
				std::vector<std::string>::iterator values = tokens.begin();
				std::vector<std::string>::iterator tokenEnd = tokens.end();
				while (++values != tokenEnd) {
					dp.push_back(atof(values->c_str()));
					cerr<<"DP: "<<*values<<"\t"<<dp[dp.size()-1]<<"\n";
				}
				cerr<< dp.size() << "Total DP values to be used.";
			}
			else {
				if (tokens[0][0] != '#') {
					std::stringstream ss(line);
					std::string pop, popFilename, word;
					ss>>pop>>popFilename;

					std::stringstream desc;
					while (!ss.eof()) {
						ss>>word;
						desc<<word;
					}
					
					splines.push_back(PopulationSpline(pop, desc.str(), popFilename));
				}
			}
		}
	}
}




int LdSplineImporter::GetPopID(soci::session& sociDB, const char* popName, float threshold) {
	int popID;
	sociDB<<"SELECT population_id FROM populations WHERE population_label=:type", soci::into(popID), soci::use(std::string(popName));
	return popID;
}


void LdSplineImporter::Process(soci::session& sociDB) {

	std::vector<PopulationSpline>::iterator spItr = splines.begin();
	std::vector<PopulationSpline>::iterator spEnd = splines.end();


	while (spItr != spEnd) {
		std::map<std::string, int> popIDs;
		InitPopulationIDs(sociDB, popIDs, *spItr, "DP", dp);
		InitPopulationIDs(sociDB, popIDs, *spItr, "RS", rs);


		LdSpline ldspline;
		ldspline.OpenBinary(spItr->filename.c_str());

		std::map<std::string, LocusLookup> chromosomes = ldspline.GetChromosomes();
		std::map<std::string, LocusLookup>::iterator chr = chromosomes.begin();
		std::map<std::string, LocusLookup>::iterator end = chromosomes.end();

		while (chr != end) {
			LoadGenes(sociDB, chr->second.Chromosome().c_str());
			ProcessLD(sociDB, chr->second, *spItr, popIDs);
			chr->second.Release();
			chr++;
		}
		spItr++;
	}
}

void LdSplineImporter::ProcessLD(soci::session& sociDB, LocusLookup& chr, PopulationSpline& sp, std::map<std::string, int>& popIDs) {
	std::vector<RegionBoundary>::iterator regItr = regions.begin();
	std::vector<RegionBoundary>::iterator regEnd = regions.end();

	cerr<<chr.Chromosome()<<"(";cerr.flush();
	std::map<std::string, int>::iterator pi = popIDs.begin();
	std::map<std::string, int>::iterator pe = popIDs.end();
	while (pi != pe) {
		cerr<<pi->first<<" ";cerr.flush();
		pi++;
	}
	int incCount = 0;
	while (regItr != regEnd) {
		int lower = regItr->lower, upper=regItr->upper;
		//cerr<<"\t--"<<regItr->geneID<<"\n";
		std::vector<float>::iterator vItr = dp.begin();
		std::vector<float>::iterator vEnd = dp.end();

		while (vItr != vEnd) {
			int popID = popIDs[sp.GetPopulationName("DP", *vItr)];

			pair<int, int> bounds = chr.GetRangeBoundariesDP(lower, upper, *vItr);
			if (bounds.first != lower || bounds.second != upper) {
				incCount++;
//				cerr<<"Update: "<<regItr->geneID<<" "<<lower<<" "<<upper<<" : "<<bounds.first<<" "<<bounds.second<<"\n";
				sociDB << "UPDATE region_bounds SET start=:start, end=:end WHERE gene_id=:geneID AND population_id=:popID", use(bounds.first), use(bounds.second), use(regItr->geneID), use(popID);
			}
			vItr++;
		}

		vItr = rs.begin();
		vEnd = rs.end();
		while (vItr != vEnd) {
			int popID = popIDs[sp.GetPopulationName("RS", *vItr)];
			pair<int, int> bounds = chr.GetRangeBoundariesRS(lower, upper, *vItr);
			if (bounds.first != lower || bounds.second != upper) {
				incCount++;
//				cerr<<"Update: "<<regItr->geneID<<" "<<lower<<" "<<upper<<" : "<<bounds.first<<" "<<bounds.second<<"\n";
				sociDB << "UPDATE region_bounds SET start=:start, end=:end WHERE gene_id=:geneID AND population_id=:popID", use(bounds.first), use(bounds.second), use(regItr->geneID), use(popID);
			}
			vItr++;
		}
		
		regItr++;

	}
	cerr<<")\t"<<incCount<<"\n";
}

void LdSplineImporter::InitPopulationIDs(soci::session& sociDB, std::map<std::string, int>& popIDs, PopulationSpline& sp, const char *type, std::vector<float>& stats) {
	std::vector<float>::iterator sItr = stats.begin();
	std::vector<float>::iterator sEnd = stats.end();


	while (sItr != sEnd) {
		cerr<<"\t"<<*sItr++<<"\n";
	}
	sItr = stats.begin();
	while (sItr != sEnd) {
		cerr<<"Initializing Population: "<<type<<", "<<*sItr<<"\t"<<stats.size()<<"\n";
		std::string popName = sp.GetPopulationName(type, *sItr);
		int popID = -1;
		sociDB<<"SELECT population_id FROM populations WHERE population_label=:name", into(popID), use(popName);
		if (popID > 0){
			cerr<<"Clearing out all bounds associated with population"<<popID<<" ("<<popName<<")\n";
			sociDB<<"DELETE FROM region_bounds WHERE population_id=:id", use(popID);
			sociDB<<"DELETE FROM populations WHERE population_id=:id", use(popID);
		} else {
			sociDB<<"SELECT MAX(population_id) FROM populations", into(popID);
			popID++;
		}
		std::stringstream itemdesc;
		itemdesc<<sp.desc<<" with "<<type<<" cutoff "<<*sItr;
		sociDB<<"INSERT INTO populations VALUES (:id, :name, :desc, :itemdesc)", use(popID), use(popName), use(sp.desc), use(itemdesc.str());
		stringstream sql;
		sql << "INSERT INTO region_bounds SELECT gene_id, "<<popID<<", start, end FROM region_bounds WHERE population_id=0";
		sociDB<<sql.str();
		popIDs[popName] = popID;
		sItr++;
	}
}


void LdSplineImporter::LoadGenes(soci::session& sociDB, const char *chrom) {
	cerr<<"SQL: "<<"SELECT gene_id, chrom, start, end FROM regions NATURAL JOIN region_bounds WHERE population_id=0 and chrom='"<<chrom<<"' ORDER BY chrom, start\n";
	rowset<row> genes = (sociDB.prepare<<"SELECT gene_id, chrom, start, end FROM regions NATURAL JOIN region_bounds WHERE population_id=0 and chrom=:chrom ORDER BY chrom, start", use(std::string(chrom)));
	regions.clear();

	for (rowset<row>::const_iterator itr = genes.begin(); itr != genes.end(); ++itr) {
		int geneID, start, end;
		std::string chrom;
		row const& row = *itr;
		geneID = row.get<int>(0);
		chrom = row.get<string>(1);
		start = row.get<int>(2);
		end   = row.get<int>(3);
		regions.push_back(RegionBoundary(geneID, chrom, start, end));
		if (start == 0) {
			cerr<<"LoadGenes: "<<chrom<<" "<<geneID<<" "<<start<<" "<<end<<"\n";
			cerr<<"--"<<regions[regions.size() - 1].lower<<"\n";
		}
	}
	cerr<<"Total Regions: "<<regions.size()<<"\n";
}

}
