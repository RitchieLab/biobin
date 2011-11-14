/* 
 * File:   Application.cpp
 * Author: torstees
 * 
 * Created on March 28, 2011, 1:45 PM
 */
#include <string>
#include <sstream>
#include "application.h"
#include <soci-sqlite3.h>
#include "utility/filetools.h"
#include "taskgenecoverage.h"
#include <iomanip>
#include "ldsplineimporter.h"
#include "liftover/converterdb.h"

namespace Biofilter {

bool Application::errorExit = false;

std::string Application::GetReportLog() {
	return reportLog.str();
}


void Application::SetVariationFilename(const char *filename) {
	sociDB<<"UPDATE versions SET version=:var WHERE element='variations'", soci::use(std::string(filename));
}


void Application::ListGenes(std::ostream& os, Utility::StringArray& aliasList, Utility::StringArray& aliasType) {
	std::string aliasClause = "";
	std::string typeClause = "";
	if (aliasList.size() > 0)
		aliasClause = std::string("( alias LIKE '") + Utility::Join(aliasList, "%' OR alias LIKE '") + "%' ) ";
	if (aliasType.size() > 0)
		typeClause = std::string("( region_alias_type_desc LIKE '") + Utility::Join(aliasType, "%' OR region_alias_type_desc LIKE '") + "%' )";
	std::string whereClause = aliasClause + typeClause;
	if (aliasClause.length() > 0 && typeClause.length() > 0)
		whereClause = std::string("AND ")  + aliasClause + " AND " + typeClause;
	else if (whereClause.length() > 0)
		whereClause = std::string("AND ") + whereClause;
	os<<"Name\tAlias\tChrom\tStart\tEnd\tDescription\tAlias Type\n";
	soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT primary_name, alias, chrom, start, end, description, region_alias_type_desc FROM regions NATURAL JOIN region_bounds NATURAL JOIN region_alias NATURAL JOIN region_alias_type WHERE population_id=0 "<<whereClause);
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		os<<row.get<std::string>(0)<<"\t"
			 <<row.get<std::string>(1)<<"\t"
			 <<row.get<std::string>(2)<<"\t"
			 <<row.get<int>(3)<<"\t"
			 <<row.get<int>(4)<<"\t"
			 <<row.get<std::string>(5)<<"\t"
			 <<row.get<std::string>(6)<<"\n";
	}
}


void Application::ListGroupIDs(std::ostream& os, Utility::StringArray& searchList) {
	std::string clause = "";

	if (searchList.size() > 0) 
		clause = std::string("WHERE ( group_name LIKE '%") + Utility::Join(searchList, "%' OR group_name LIKE '%") + "%' OR group_desc LIKE '%" + Utility::Join(searchList, "%' OR group_desc LIKE '%'") + "%') ";
	os<<"ID\tName\tDescription\n";
	soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT group_id, group_name, group_desc FROM groups "<<clause);
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		os<<row.get<int>(0)<<"\t"
			 <<row.get<std::string>(1)<<"\t"
			 <<row.get<std::string>(2)<<"\n";
	}

}

void Application::ListPopulationIDs(std::ostream& os) {
	soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT population_label, pop_description FROM populations");
	os<<"Label\tDescription\n";
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		os<<row.get<std::string>(0)<<"\t"<<row.get<std::string>(1)<<"\n";
	}
}

//Gene coverage is a special report that runs independent of all the other reports
//So, it is basically a completely separate path from the rest. So, we will be using
//stuff without regard to other possible executions-so don't forget that
void Application::GeneCoverage(Utility::StringArray &rsSources, Utility::StringArray &mapSources, const char *geneFilename, const char *population) {
	Utility::StringArray missingAliases;
	Utility::StringArray aliasList;
	
	if (strcmp(geneFilename, "ALL") != 0) {
		std::string fileContents = Utility::LoadContents(geneFilename);
		aliasList = Utility::Split(fileContents.c_str());
	}


	Task::GeneCoverage gc(&dataset, &regions, variationFilename.c_str());
	gc.AddSources(rsSources, mapSources, &buildConverter);
	//Selectively load genes
	LoadRegionData(population, missingAliases, aliasList);
	std::string reportFilename = AddReport("gene-coverage", "csv", "Gene Coverage Report");
	gc.GenerateTxtReport(reportFilename.c_str());
	
}


std::multimap<uint, uint> Application::BuildSnpGeneMap() {
	regions.BuildSnpGeneMap(geneLookup);
	return geneLookup;
}

std::string Application::AddReport(const char *suffix, const char *extension, const char *description) {
	std::string sfx = "";
	if (strlen(suffix) > 0)
		sfx = "-"+std::string(suffix);
	std::string newFilename = reportPrefix + sfx + "." + extension;
	reportLog<<std::setw(50)<<std::right<<newFilename<<" : "<<description<<"\n";
	return newFilename;
}

void Application::ProduceModels(std::ostream& os) {
	std::map<uint, Knowledge::GroupManagerDB>::iterator itr = groups.begin();
	std::map<uint, Knowledge::GroupManagerDB>::iterator end = groups.end();

	while (itr != end) 
		itr++->second.GenerateGeneGeneModels(geneGeneModels, regions, os);
}

void Application::LoadBuildConverter(const char *build) {
	buildConverter.LoadFromDB(build, sociDB);
}

uint Application::LoadMapData(const char *filename, const char *genomicBuild, Knowledge::SnpDataset& lostSNPs, bool performAlignment) {
	std::ifstream file(filename);
	if (!file.good())
		throw Utility::Exception::FileNotFound(filename);

	char line[4096];
	LiftOver::ConverterDB cnv;
	int chainCount = cnv.LoadFromDB(genomicBuild, sociDB);

	if (chainCount>=0) {
		LiftOver::SnpArray snps;

		while (file.good()) {
			file.getline(line, 4096);
			std::stringstream ss(line);
			Utility::StringArray words = Utility::Split(line);
			
			if (words.size() > 3)
				snps.push_back(LiftOver::SNP(Utility::ChromToInt(words[0].c_str()), atoi(words[3].c_str()), words[1].c_str()));
				//dataset.AddSNP(Utility::ChromToInt(words[0].c_str()), atoi(words[3].c_str()), words[1].c_str());
		}

		std::multimap<LiftOver::SNP, LiftOver::SNP> converted;
		cnv.ConvertDataset(snps, converted);
		std::multimap<LiftOver::SNP, LiftOver::SNP>::iterator itr = converted.begin();
		std::multimap<LiftOver::SNP, LiftOver::SNP>::iterator end = converted.end();
		
		std::stringstream missingSNPs;
		while (itr != end) {
			//Report any bad SNPS
			if (itr->second.pos == 0) 
				missingSNPs<<itr->first.RSID()<<"\t"<<Utility::ChromFromInt(itr->first.chrom)<<"\t"<<itr->first.pos<<"\n";
			else 
				dataset.AddSNP(itr->second.chrom, itr->second.pos, itr->second.RSID().c_str());
			itr++;
		}
		
		if (missingSNPs.str().length() > 0) {
			std::cerr<<"SNPs that weren't translated properly to the new build:\n";
			std::cerr<<"Chr/tPos/tRS ID\n"<<missingSNPs.str()<<"\n";
		}
		
		if (performAlignment)
			dataset.AlignData();
		return dataset.Size();
	}
	return 0;
}

uint Application::LoadSnpsSource(const char *filename, std::set<std::string>& lostSNPs) {
	std::string fileContents = Utility::LoadContents(filename);
	std::set<std::string> rsIDs;
	Utility::StringArray lines = Utility::Split(fileContents.c_str(), "\n");
	Utility::StringArray::iterator itr = lines.begin();
	Utility::StringArray::iterator end = lines.end();
	while (itr != end) {
		std::string rsid = Utility::Locus::_RSID(itr++->c_str());
		rsIDs.insert(rsid);
	}
	
	std::set<std::string> snpsFound;
	uint count = dataset.LoadData(rsIDs, snpsFound);

	std::set_difference(rsIDs.begin(), rsIDs.end(), snpsFound.begin(), snpsFound.end(), std::inserter(lostSNPs, lostSNPs.begin()));
	
	return count;
}


uint Application::GetPopulationID(const char *pop) {
	int popID = 0;

	sociDB<<"SELECT population_id FROM populations WHERE population_label='"<<pop<<"'", soci::into(popID);
	return popID;
}

void Application::LoadLdSpline(const char *cfg) {
	LdSplineImporter splineMgr;
	splineMgr.LoadConfiguration(cfg);
	splineMgr.Process(sociDB);
}

uint Application::LoadRegionData(const char *pop, Utility::StringArray& aliasesNotFound, Utility::StringArray& aliasList) {
	std::map<std::string, uint> idLookup;
	uint popID = GetPopulationID(pop);

	regions.LoadFromDB(sociDB, popID, aliasList, idLookup);
	regions.AssociateSNPs(dataset);

	Utility::StringArray::iterator itr = aliasList.begin();
	Utility::StringArray::iterator end = aliasList.end();
	while (itr != end) {
		if (idLookup.find(*itr) == idLookup.end())
			aliasesNotFound.push_back(*itr);
		itr++;
	}
			
	return idLookup.size();
}


/**
 * This function depends on regions having been properly loaded
 */
uint Application::LoadGroupDataByName(Utility::StringArray& userDefinedGroups, 
		Utility::StringArray& groupNames,
		Utility::IdCollection& groupIDs) {
	int maxGroupID;
	sociDB<<"SELECT MAX(group_id) FROM groups", soci::into(maxGroupID);
	
	int groupTypeCount;
	sociDB<<std::string("SELECT COUNT(group_type_id) FROM group_type WHERE group_type_id IN (")+ Utility::Join(groupIDs, ",") + ") ", soci::into(groupTypeCount);
	uint idCount = groupIDs.size();
	std::string sql = "SELECT group_type_id, group_type, role_id FROM group_type";
	if (groupTypeCount > 0) {
		sql = std::string("SELECT group_type_id, group_type, role_id FROM group_type WHERE group_type_id IN (") + Utility::Join(groupIDs, ",") + ") ";
	}
	soci::rowset<soci::row> rs = (sociDB.prepare << sql);

	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		uint id = row.get<int>(0);
		std::string name = row.get<std::string>(1);
		Knowledge::MetaGroup::Type roleID = row.get<int>(2);
		groupIDs.erase(id);
		groups[id] = Knowledge::GroupManagerDB(id, roleID, name.c_str(), "");
	}

	Utility::StringArray unmatchedAliases;
	Utility::StringArray::iterator udItr = userDefinedGroups.begin();
	Utility::StringArray::iterator udEnd = userDefinedGroups.end();

	uint totalGroupsLoaded = 0;
	while (udItr != udEnd) {
		//Give some bogus groupType, since it will be found in the file
		Knowledge::GroupManagerDB udGroup(++maxGroupID, 1, udItr->c_str());
		std::map<std::string, uint> nameToId = udGroup.LoadArchive(udItr->c_str(), regions, maxGroupID + 1, unmatchedAliases);
		groups[udGroup.id] = udGroup;
		totalGroupsLoaded += nameToId.size();
		maxGroupID += nameToId.size();
		udItr++;
	}

	std::map<uint, Knowledge::GroupManagerDB>::iterator itr = groups.begin();
	std::map<uint, Knowledge::GroupManagerDB>::iterator end = groups.end();

	while (itr != end) {
		uint count = itr->second.LoadFromDB(sociDB, groupIDs, regions, groupNames);
		if (count > 0) {
			std::cerr<<itr++->second.name<<"\t"<<count<<"\n";
			totalGroupsLoaded += count;
		} else
			groups.erase(itr++);
	}

	return totalGroupsLoaded;

}



void Application::InitBiofilter(const char *filename, bool reportVersion) {
	dbFilename = filename;
	boost::filesystem::path dbPath = boost::filesystem::path(dbFilename);
	bool fileFound = false;
	if (boost::filesystem::is_regular_file(dbPath)) {
		fileFound = true;
	}else{
		#ifdef DATA_DIR
			if (dbPath.is_relative()){
				dbPath = (boost::filesystem::path(std::string(DATA_DIR))/=(dbPath));
				if (boost::filesystem::is_regular_file(dbPath)){
					fileFound=true;
				}
			}
		#endif
	}

	if (!fileFound){
		throw Utility::Exception::FileNotFound(filename);
	}

	try {
		std::string cnxParam = "dbname="+dbPath.native()+" timeout=2500";
		sociDB.open(soci::sqlite3, cnxParam.c_str());
		std::string dbSnp, ensembl, hapmap, build, variations;
		sociDB<<"SELECT version FROM versions WHERE element='ncbi'", soci::into(dbSnp);
		sociDB<<"SELECT version FROM versions WHERE element='ensembl'", soci::into(ensembl);
		sociDB<<"SELECT version FROM versions WHERE element='hapmap'", soci::into(hapmap);
		sociDB<<"SELECT version FROM versions WHERE element='variations'", soci::into(variationFilename);
		sociDB<<"SELECT version FROM versions WHERE element='build'", soci::into(build);
		this->varVersion				= atoi(variations.c_str());

		if (reportVersion) {
			std::cerr<<"\n------------------------- Dependency Versions ----------\n";
			std::cerr<<std::setw(38)<<std::right<<"dbSNP : "<<dbSnp<<"\n";
			std::cerr<<std::setw(38)<<std::right<<"Ensembl : "<<ensembl<<"\n";
			std::cerr<<std::setw(38)<<std::right<<"Hap Map LD : "<<hapmap<<"\n";
			std::cerr<<std::setw(38)<<std::right<<"Variation Filename : "<<variationFilename<<"\n";
			std::cerr<<std::setw(38)<<std::right<<"Genome Build : "<<build<<"\n";

		}
		dataset.SetVariationsFilename(variationFilename.c_str());

	} catch (soci::soci_error const &e) {
		std::cerr<<"Problems were encountered trying to open the database, "<<dbFilename<<". Error: "<<e.what()<<"\n";
	}


	try {
		soci::rowset<soci::row> rs = (sociDB.prepare << "SELECT id, role FROM snp_role");

		for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			soci::row const& row = *itr;
			uint id = row.get<int>(0);
			std::string name = row.get<std::string>(1);

			dataset.RoleDescription(id, name.c_str());
		}
	} catch (soci::soci_error const &e) {
		std::cerr<<"An error was encountered trying to load the role data. SNP role information will not be available.\n";
	}
}

}
