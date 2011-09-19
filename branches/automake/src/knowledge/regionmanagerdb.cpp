/* 
 * File:   dbregionmanager.cpp
 * Author: torstees
 * 
 * Created on March 9, 2011, 2:45 PM
 */

#include "regionmanagerdb.h"
#include <soci-sqlite3.h>
#include "utility/locus.h"

namespace Knowledge {


void RegionManagerDB::GenerateLookupTable(std::map<uint, uint>& lookup) {
	uint count = regions.size();
	for (uint i=0; i<count; i++)
		lookup[regions[i].id] = i;
}






void RegionManagerDB::LoadRegionAliases(soci::session& sociDB, Utility::StringArray& aliasList, std::map<std::string, uint>& aliasToID) {
	std::string conditions = "";
	if (aliasList.size() > 0)
		conditions = " WHERE alias IN (" + Utility::Join(aliasList, ",") + ")";

	soci::rowset<soci::row> rs = (sociDB.prepare <<"SELECT gene_id, alias FROM region_alias"<<conditions);

	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		uint geneID			= row.get<int>(0);
		std::string alias	= row.get<std::string>(1);
		aliasToID[alias] = geneID;
	}
}


uint RegionManagerDB::LoadFromDB(soci::session& sociDB, uint popID) {
	Utility::IdCollection emptyIDs;
	return LoadFromDB(sociDB, popID, emptyIDs);
}

uint RegionManagerDB::LoadFromDB(soci::session& sociDB, uint popID, Utility::StringArray& aliasList, std::map<std::string, uint>& aliasToID) {
	LoadRegionAliases(sociDB, aliasList, aliasToID);
	std::map<std::string, uint>::iterator itr = aliasToID.begin();
	std::map<std::string, uint>::iterator end = aliasToID.end();

	Utility::IdCollection ids;
	if (aliasList.size() > 0) {
		while (itr != end)
			ids.insert(itr++->second);
	}
	return LoadFromDB(sociDB, popID, ids);
}

uint RegionManagerDB::LoadFromDB(soci::session& sociDB, uint popID, const Utility::IdCollection& ids) {
	std::string fromIDs = "";
	if (ids.size() > 0) {
		std::string idlist = Utility::Join<std::set<uint> >(ids, ",");
		fromIDs = std::string("AND regions.gene_id IN (") + idlist + std::string(")");
	}

	int geneID, eStart, eStop, oStart, oStop;
	soci::indicator aliasIndicator;
	std::string name, aliases, chr;
	//Something about joining region_bounds to region_bounds (with different population IDs) was killing the query time, so
	//I'm splitting these two up. The query ran just fine in sqlite manager, so I have no idea why it took orders of magnitude
	//longer through soci. 
	soci::statement st = (sociDB.prepare <<"SELECT regions.gene_id, primary_name, chrom, eff.start, eff.end, group_concat(alias) "<<
		 "FROM regions NATURAL JOIN region_bounds eff LEFT JOIN region_alias ON (regions.gene_id=region_alias.gene_id) "<<
		 "WHERE eff.population_id="
		 <<popID<<" "<<fromIDs<<"  GROUP BY regions.gene_id", soci::into(geneID), soci::into(name), soci::into(chr), soci::into(eStart), soci::into(eStop), soci::into(aliases, aliasIndicator));

	st.execute();
	while (st.fetch()) {
		char chrom			= Utility::ChromToInt(chr.c_str());
		std::string aliasList = "";

		if (aliasIndicator != soci::i_null)
			aliasList = aliases.c_str();

		Region &r = AddRegion(name.c_str(), geneID, Utility::ChromToInt(chr.c_str()), eStart, eStop, aliasList.c_str());
		r.chrom = chrom;
	}

	st = (sociDB.prepare<<"SELECT gene_id, start, end FROM region_bounds WHERE population_id=0", soci::into(geneID), soci::into(oStart), soci::into(oStop));
	st.execute();

	uint missing = (uint)-1;
	while (st.fetch()) {
		uint id = (*this)(geneID);
		if (id != missing) {
			regions[id].trueStart = oStart;
			regions[id].trueEnd   = oStop;
		}
	}
	std::cerr<<regions.size()<<" "<<"Regions Loaded: "<<"\n";

	return regions.size();
}

/*
uint RegionManagerDB::LoadFromDB(soci::session& sociDB, uint popID, const Utility::IdCollection& ids) {
	std::string fromIDs = "";
	if (ids.size() > 0) {
		std::string idlist = Utility::Join<std::set<uint> >(ids, ",");
		fromIDs = std::string("AND regions.gene_id IN (") + idlist + std::string(")");
	}
	std::cout<<"SELECT regions.gene_id, primary_name, chrom, eff.start, eff.end, orig.start, orig.end, group_concat(alias) "<<
		 "FROM region_bounds orig NATURAL JOIN regions NATURAL JOIN region_bounds eff LEFT JOIN region_alias ON (regions.gene_id=region_alias.gene_id) "<<
		 "WHERE orig.population_id=0 AND eff.population_id="
		 <<popID<<" "<<fromIDs<<"  GROUP BY regions.gene_id"<<"\n";
	soci::rowset<soci::row> rs = (sociDB.prepare <<"SELECT regions.gene_id, primary_name, chrom, eff.start, eff.end, orig.start, orig.end, '' "
		 "FROM region_bounds orig NATURAL JOIN regions NATURAL JOIN region_bounds eff  "
		 "WHERE orig.population_id=0 AND eff.population_id="
		 <<popID<<" "<<fromIDs<<"  GROUP BY regions.gene_id");
	//soci::rowset<soci::row> rs = (sociDB.prepare <<"SELECT regions.gene_id, primary_name, chrom, eff.start, eff.end, orig.start, orig.end, group_concat(alias) FROM region_bounds orig NATURAL JOIN regions NATURAL JOIN region_bounds eff NATURAL JOIN region_alias WHERE orig.population_id=0 AND eff.population_id="<<popID<<" "<<fromIDs);
	regions.clear();
	for (soci::rowset<soci::row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		soci::row const& row = *itr;
		uint geneID			= row.get<int>(0);
		std::string name	= row.get<std::string>(1);
		char chrom			= Utility::ChromToInt(row.get<std::string>(2).c_str());
		uint effstart		= row.get<int>(3);
		uint effstop		= row.get<int>(4);
		uint origstart		= row.get<int>(5);
		uint origstop		= row.get<int>(6);
		std::string aliases = "";		//row.get<std::string>(7, "");
		Region &r = AddRegion(name.c_str(), geneID, effstart, effstop, origstart, origstop, aliases.c_str());
		r.chrom = chrom;
	}

	std::cerr<<regions.size()<<" "<<"Regions Loaded: "<<"\n";

	return regions.size();

}
*/
}



#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

class DbRegionManagerTest : public ::testing::Test {
public:
	Knowledge::RegionManagerDB regions;
	SnpDataset dataset;


	virtual void SetUp() {
		//OK, let's create the database
		remove("fake-bio.db");
		soci::session sociDB;
		sociDB.open(soci::sqlite3, "dbname=fake-bio.db timeout=2500");
		sociDB<<"CREATE TABLE regions (gene_id INTEGER UNIQUE PRIMARY KEY, primary_name VARCHAR(128) UNIQUE, chrom VARCHAR(2), description TEXT)";
		sociDB<<"CREATE TABLE region_bounds (gene_id INTEGER, population_id INTEGER, start INTEGER, end INTEGER)";
		sociDB<<"CREATE TABLE region_alias_type (region_alias_type_id INTEGER UNIQUE PRIMARY KEY, region_alias_type_desc TEXT)";
		sociDB<<"CREATE TABLE region_alias (region_alias_type_id INTEGER, alias VARCHAR(64), gene_id INTEGER, gene_count INTEGER, UNIQUE(gene_id, alias, region_alias_type_id))";
		sociDB<<"INSERT INTO regions VALUES (1, 'A1BG', '19', 'alpha-1-B glycoprotein')";
		sociDB<<"INSERT INTO regions VALUES (2, 'A2M',  '12', 'alpha-2-macroglobulin')";
		sociDB<<"INSERT INTO regions VALUES (3, 'A2MP1', '12', 'aplha-2-macroglobulin pseudogene')";
		sociDB<<"INSERT INTO regions VALUES (9,'NAT1','8','N-acetyltransferase 1 (arylamine N-acetyltransferase)')";
		sociDB<<"INSERT INTO regions VALUES (10,'NAT2','8','N-acetyltransferase 2 (arylamine N-acetyltransferase)')";
		sociDB<<"INSERT INTO regions VALUES (11, 'AACCP', '8', 'arylamide acetylase pseudogene')";
		sociDB<<"INSERT INTO region_alias VALUES (2000, 'A8K052', 1, 1)";
		sociDB<<"INSERT INTO region_alias VALUES (1300, 'ABG', 1, 1)";
		sociDB<<"INSERT INTO region_alias VALUES (2000, 'C9J773', 2, 1)";
		sociDB<<"INSERT INTO region_alias VALUES (1300, 'CPAMD5', 2, 1)";

		sociDB<<"INSERT INTO region_bounds VALUES ('1','0','58858171','58864864')";
		sociDB<<"INSERT INTO region_bounds VALUES ('2','0','9220303','9268557')";
		sociDB<<"INSERT INTO region_bounds VALUES ('3','0','9384201','9386907')";
		sociDB<<"INSERT INTO region_bounds VALUES ('9','0','18067289','18081197')";
		sociDB<<"INSERT INTO region_bounds VALUES ('10','0','18248754','18258722')";
		sociDB<<"INSERT INTO region_bounds VALUES ('11','0','18227418','18229386')";
		sociDB<<"INSERT INTO region_bounds VALUES ('1','1','58857171','58866864')";
		sociDB<<"INSERT INTO region_bounds VALUES ('2','1','9220003','9268657')";
		sociDB<<"INSERT INTO region_bounds VALUES ('3','1','9382201','9387907')";
		sociDB<<"INSERT INTO region_bounds VALUES ('9','1','18067289','18081197')";
		sociDB<<"INSERT INTO region_bounds VALUES ('10','1','18228754','18258762')";		// Overlaps with gene 11
		sociDB<<"INSERT INTO region_bounds VALUES ('11','1','18027418','18229586')";



		//Add 3 SNPs to chromosome 19 all interior to A1BG
		dataset.AddSNP(19, 58858171, "rs1", 1);	// Interior
		dataset.AddSNP(19, 58860864, "rs2", 1);	// Interio
		dataset.AddSNP(19, 58864864, "rs3", 1);	// Interior

		//Add 5 SNPs to 12 - 4 in A2M (only one interior), 0 in A2MP1
		dataset.AddSNP(12, 9220305, "rs4", 1);		// Interior
		dataset.AddSNP(12, 9220603, "rs5", 1);		// LD
		dataset.AddSNP(12, 9268657, "rs6", 1);		// LD
		dataset.AddSNP(12, 9266657, "rs7", 1);		// LD
		dataset.AddSNP(12, 9000000, "rs8", 1);		// not in gene

		//Add 10 SNPs with 2 in overlapping region
		//4 in NAT2 only (2 interior)
		dataset.AddSNP(8, 18256754, "rs9", 1);		// Interior
		dataset.AddSNP(8, 18258722, "rs10", 1);	// Interior Edge
		dataset.AddSNP(8, 18239386, "rs11", 1);	// LD
		dataset.AddSNP(8, 18258762, "rs12", 1);	// LD Edge

		//3 Overlapping
		dataset.AddSNP(8, 18228754, "rs13", 1);	// LD edge (NAT2)
		dataset.AddSNP(8, 18229386, "rs14", 1);	// LD Edge (AACCP)
		dataset.AddSNP(8, 18229586, "rs16", 1);	// LD edge

		//3 in AACCP (2 interior)
		dataset.AddSNP(8, 18227418, "rs15", 1);	// interior edge
		dataset.AddSNP(8, 18027418, "rs17", 1);	// LD Edge
		dataset.AddSNP(8, 18027618, "rs18", 1);   // Interior


	}

	virtual void TearDown() {
		remove("fake-bio.db");
	}

};

TEST_F(DbRegionManagerTest, AliasTest) {
	soci::session sociDB;
	sociDB.open(soci::sqlite3, "dbname=fake-bio.db timeout=2500");
	Knowledge::RegionManagerDB regions;
	regions.LoadFromDB(sociDB, 1);
	regions.AssociateSNPs(dataset);

	EXPECT_EQ("A1BG", regions["A8K052"].name);
	EXPECT_EQ("A1BG", regions["ABG"].name);
	EXPECT_EQ((uint)1, regions["A8K052"].id);
	EXPECT_EQ(19, regions["A8K052"].chrom);
	EXPECT_TRUE(regions["A8K052"].IsPresent(0));
	EXPECT_TRUE(regions["A8K052"].IsPresent(1));
	EXPECT_TRUE(regions["A8K052"].IsPresent(2));

	EXPECT_EQ("A2M", regions["C9J773"].name);
	EXPECT_EQ("A2M", regions["CPAMD5"].name);
	EXPECT_EQ(2, regions["C9J773"].id);
	EXPECT_EQ(12, regions["C9J773"].chrom);
	EXPECT_TRUE(regions["C9J773"].IsPresent(3));
	EXPECT_TRUE(regions["C9J773"].IsPresent(4));
	EXPECT_TRUE(regions["C9J773"].IsPresent(5));
}

TEST_F(DbRegionManagerTest, LookupTable) {
	soci::session sociDB;
	sociDB.open(soci::sqlite3, "dbname=fake-bio.db timeout=2500");
	Knowledge::RegionManagerDB regions;
	regions.LoadFromDB(sociDB, 1);
	regions.AssociateSNPs(dataset);

	std::map<uint, uint> lookup;
	regions.GenerateLookupTable(lookup);
	//		0=>1, 1=>2, 2=>10, 3=>11
	EXPECT_EQ(0, lookup[1]);
	EXPECT_EQ(1, lookup[2]);
	EXPECT_EQ(2, lookup[10]);
	EXPECT_EQ(3, lookup[11]);
}

TEST_F(DbRegionManagerTest, SqueezeTest) {
	soci::session sociDB;
	sociDB.open(soci::sqlite3, "dbname=fake-bio.db timeout=2500");
	Knowledge::RegionManagerDB regions;
	uint count = regions.LoadFromDB(sociDB, 1);
	EXPECT_EQ(6, count);
	EXPECT_EQ(6, regions.Size());
	
	EXPECT_EQ("A1BG", regions[0].name);
	EXPECT_EQ("A2M", regions[1].name);
	EXPECT_EQ("A2MP1", regions[2].name);
	EXPECT_EQ("NAT1", regions[3].name);
	EXPECT_EQ("NAT2", regions[4].name);
	EXPECT_EQ("AACCP", regions[5].name);

	EXPECT_EQ(1, regions[0].id);
	EXPECT_EQ(2, regions[1].id);
	EXPECT_EQ(11, regions[5].id);

	EXPECT_EQ(19, regions[0].chrom);
	EXPECT_EQ(12, regions[1].chrom);
	EXPECT_EQ(8, regions[5].chrom);

	EXPECT_EQ(58858171, regions[0].trueStart);
	EXPECT_EQ(58864864, regions[0].trueEnd);
	EXPECT_EQ(58857171, regions[0].effStart);
	EXPECT_EQ(58866864, regions[0].effEnd);

	EXPECT_EQ(18227418, regions[5].trueStart);
	EXPECT_EQ(18229386, regions[5].trueEnd);
	EXPECT_EQ(18027418, regions[5].effStart);
	EXPECT_EQ(18229586, regions[5].effEnd);

	regions.AssociateSNPs(dataset);

	EXPECT_EQ(4, regions.Size());
	EXPECT_EQ(3, regions[0].SnpCount());
	EXPECT_EQ(4, regions[1].SnpCount());
	EXPECT_EQ(7, regions[2].SnpCount());
	EXPECT_EQ(6, regions[3].SnpCount());

	EXPECT_TRUE(regions[0].IsPresent(0));
	EXPECT_TRUE(regions[0].IsPresent(1));
	EXPECT_TRUE(regions[0].IsPresent(2));

	EXPECT_TRUE(regions[1].IsPresent(3));
	EXPECT_TRUE(regions[1].IsPresent(4));
	EXPECT_TRUE(regions[1].IsPresent(5));
	EXPECT_TRUE(regions[1].IsPresent(6));

	EXPECT_TRUE(regions[2].IsPresent(8));
	EXPECT_TRUE(regions[2].IsPresent(9));
	EXPECT_TRUE(regions[2].IsPresent(10));
	EXPECT_TRUE(regions[2].IsPresent(11));
	EXPECT_TRUE(regions[2].IsPresent(12));
	EXPECT_TRUE(regions[2].IsPresent(13));
	EXPECT_TRUE(regions[2].IsPresent(14));

	EXPECT_TRUE(regions[3].IsPresent(12));
	EXPECT_TRUE(regions[3].IsPresent(13));
	EXPECT_TRUE(regions[3].IsPresent(14));
	EXPECT_TRUE(regions[3].IsPresent(15));
	EXPECT_TRUE(regions[3].IsPresent(16));
	EXPECT_TRUE(regions[3].IsPresent(17));
}

#endif //TEST_APP
