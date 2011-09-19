/* 
 * File:   genegenemodel.cpp
 * Author: torstees
 * 
 * Created on March 8, 2011, 10:29 AM
 */

#include "genegenemodel.h"
#include <soci-sqlite3.h>
#include "utility/locus.h"
#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

class GeneGeneModelTest : public ::testing::Test {
public:
	GeneGeneModel m1,m2,m3,m4,dead;
	SnpDataset dataset;
	soci::session sociDB;
	Knowledge::RegionManagerDB regions;
	
	GeneGeneModelTest() :
			m1(0, 1, 1.0),
			m2(0, 2, 1.0),
			m3(1, 3, 5.2),
			m4(2, 3, 2.0),
			dead(1, 1, 0.0) {}
	
	virtual void SetUp() {
		remove("fake-bio.db");
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

		
		// Build up the dataset
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

		regions.LoadFromDB(sociDB, 1);
		regions.AssociateSNPs(dataset);


	}

	virtual void TearDown() {
		remove("fake-bio.db");
	}
};

/**
 	GeneGeneModelTest() :
			m1(0, 1, 1.0),
			m2(0, 2, 1.0),
			m3(1, 3, 5.2),
			m4(2, 3, 2.0),
			dead(1, 1, 0.0) {}

 */
TEST_F(GeneGeneModelTest, Basics) {


	//Double check the contents of the indexes
	EXPECT_EQ((uint)2, m1.Size());
	EXPECT_EQ((uint)2, m4.Size());
	EXPECT_EQ((uint)1, dead.Size());
	EXPECT_EQ(0, m1[0]);
	EXPECT_EQ(1, m1[1]);
	EXPECT_EQ(1, m3[0]);
	EXPECT_EQ(3, m3[1]);
	EXPECT_EQ(2, m4[0]);
	EXPECT_EQ(3, m4[1]);

	float n = 0.1;
	EXPECT_EQ((float)0.1, n);

	EXPECT_EQ(1.0, m1.implicationIndex);
	EXPECT_EQ(1.0, m2.implicationIndex);
	EXPECT_EQ((float)5.2, m3.implicationIndex);
	EXPECT_EQ(2.0, m4.implicationIndex);

	//Check the cases where II is the same
	EXPECT_TRUE(m1 < m2);
	EXPECT_FALSE(m2 < m1);

	EXPECT_TRUE(m3 < m1);
	EXPECT_TRUE(m3 < m4);
	EXPECT_TRUE(m4 < m2);
	EXPECT_TRUE(m1 < m2);
	EXPECT_FALSE(m2 < m4);

	GeneGeneModel m5(m1);
	EXPECT_EQ(2, m5.Size());
	EXPECT_EQ(0, m5[0]);
	EXPECT_EQ(1, m5[1]);

	EXPECT_TRUE(m1 == m5);
	EXPECT_FALSE(m2 == m5);
}

TEST_F(GeneGeneModelTest, ModelGeneration) {
	BinaryArchive = false;
	SnpSnpModel::Collection models;
	uint count = m1.GenerateModels(models, regions);
	EXPECT_EQ(12, count);
	EXPECT_EQ(12, models.size());

	SnpSnpModel::Collection::const_iterator itr = models.begin();
	EXPECT_EQ(2, itr->Size());
	EXPECT_EQ(0, (*itr)[0]);
	EXPECT_EQ(3, (*itr)[1]);
	EXPECT_FLOAT_EQ(1.0, itr++->ImplicationIndex());
	EXPECT_EQ("0	4	1", itr++->ToString());
	EXPECT_EQ("0	5	1", itr++->ToString());
	EXPECT_EQ("0	6	1", itr++->ToString());
	EXPECT_EQ("1	3	1", itr++->ToString());
	EXPECT_EQ("1	4	1", itr++->ToString());
	EXPECT_EQ("1	5	1", itr++->ToString());
	EXPECT_EQ("1	6	1", itr++->ToString());
	EXPECT_EQ("2	3	1", itr++->ToString());
	EXPECT_EQ("2	4	1", itr++->ToString());
	EXPECT_EQ("2	5	1", itr++->ToString());
	EXPECT_EQ("2	6	1", itr->ToString());

}

TEST_F(GeneGeneModelTest, ModelEstimations) {
	GeneGeneModel m5(m1);
	EXPECT_EQ(12, m1.EstimateModelCount(regions));
	EXPECT_EQ(21, m2.EstimateModelCount(regions));
	EXPECT_EQ(24, m3.EstimateModelCount(regions));
	EXPECT_EQ(42, m4.EstimateModelCount(regions));
	EXPECT_EQ(12, m5.EstimateModelCount(regions));
}

TEST_F(GeneGeneModelTest, ModelArchive) {
	GeneGeneModel m5(m1);

	std::stringstream ss;
	m1.Write(ss, false);
	EXPECT_EQ("0	1	1\n", ss.str());
	dead.Load(ss, false);
	EXPECT_TRUE(dead == m1);
	ss.str("");

	m2.Write(ss, false);
	EXPECT_EQ("0	2	1\n", ss.str());
	dead.Load(ss, false);
	EXPECT_TRUE(dead == m2);

	ss.str("");
	m3.Write(ss, false);
	EXPECT_EQ("1	3	5.2\n", ss.str());
	dead.Load(ss, false);
	EXPECT_TRUE(dead == m3);

	ss.str("");
	m4.Write(ss, false);
	EXPECT_EQ("2	3	2\n", ss.str());
	dead.Load(ss, false);
	EXPECT_TRUE(dead == m4);

	ss.str("");
	m3.Write(ss, true);
	std::stringstream inss(ss.str());
	dead.Load(inss, true);

	EXPECT_TRUE(dead == m3);

	ss.str("");
	m5.Write(ss, false);
	EXPECT_EQ("0	1	1\n", ss.str());

}




#endif //TEST_APP
