/* 
 * File:   regionmanager.cpp
 * Author: torstees
 * 
 * Created on March 7, 2011, 2:33 PM
 */

#include "regionmanager.h"


Knowledge::ModelGenerationMode::Type Knowledge::RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::ALL_MODELS;

#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

class RegionManagerTest : public ::testing::Test {
public:
	Knowledge::RegionManager regions;


	virtual void SetUp() {

		{
		Region &r = regions.AddRegion("r1", 1, 0, 25, 10, 15);
		r.AddMetaIDs(MetaGroup::DiseaseIndependent, Utility::ToSet<uint>("1,2,3", ","));
		r.AddMetaIDs(MetaGroup::DiseaseDependent, Utility::ToSet<uint>("4,5", ","));
		r.AddAliases("reg1,1");
		r.AddSNPs(Utility::ToSet<uint>("100,115,125,200,250,255,1000", ","));
		}

		{
		Region &r = regions.AddRegion("r2", 2, 20, 50, 30, 40);
		r.AddMetaIDs(MetaGroup::DiseaseIndependent, Utility::ToSet<uint>("1,3", ","));
		//r.AddMetaIDs(MetaGroup::DiseaseDependent, Utility::ToSet<uint>("5", ","));
		//r.AddMetaIDs(MetaGroup::DiseaseDependent, "");
		r.AddAliases("region2,02,x2");
		r.AddSNPs(Utility::ToSet<uint>("300 314 325 12333 51212"));
		}

		{
		Region &r = regions.AddRegion("r3", 3, 100, 150);
		r.AddMetaIDs(MetaGroup::DiseaseIndependent, Utility::ToSet<uint>("2"));
		r.AddMetaIDs(MetaGroup::DiseaseDependent, Utility::ToSet<uint>("4,6", ","));
		r.AddAliases("Number3,three");
		r.AddSNPs(Utility::ToSet<uint>("1234 2345 345678"));
		}
	}

	virtual void TearDown() { }

};

TEST_F(RegionManagerTest, SqueezeTest) {
	regions.AddRegion("r4", 4, 5, 10);
	regions.AddRegion("r5", 5, 10, 20);
	EXPECT_EQ(5, regions.Size());
	regions.Squeeze();
	EXPECT_EQ(3, regions.Size());
	Region r = regions[0];
	EXPECT_EQ("r1", r.name);
	EXPECT_EQ(1, r.id);
	EXPECT_EQ(0, r.effStart);
	EXPECT_EQ(25, r.effEnd);
	EXPECT_EQ(10, r.trueStart);
	EXPECT_EQ(15, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ("reg1", r.aliases[0]);
	EXPECT_EQ("1", r.aliases[1]);
	EXPECT_EQ(3, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());
	EXPECT_EQ(2, r.CountDDCapable());

	r = regions[1];
	EXPECT_EQ("r2", r.name);
	EXPECT_EQ(2, r.id);
	EXPECT_EQ(20, r.effStart);
	EXPECT_EQ(50, r.effEnd);
	EXPECT_EQ(30, r.trueStart);
	EXPECT_EQ(40, r.trueEnd);
	EXPECT_EQ(0, r.CountDDCapable());

	r = regions[2];
	EXPECT_EQ("r3", r.name);
	EXPECT_EQ(3, r.id);
	EXPECT_EQ(100, r.effStart);
	EXPECT_EQ(150, r.effEnd);
	EXPECT_EQ(100, r.trueStart);
	EXPECT_EQ(150, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(1, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());
	EXPECT_EQ(2, r.CountDDCapable());

	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::DD_ONLY;
	Utility::IdCollection ids = Utility::ToSet<uint>("0,1,2", ",");
	EXPECT_TRUE(regions.DoGenerateModels(ids));
	ids = Utility::ToSet<uint>("1,1", ",");
	EXPECT_FALSE(regions.DoGenerateModels(ids));

	EXPECT_TRUE(regions.ValidGeneGene(0,1));
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::GROUP_LEVEL;
	EXPECT_FALSE(regions.DoGenerateModels(ids));
	ids = Utility::ToSet<uint>("1,2", ",");
	EXPECT_TRUE(regions.DoGenerateModels(ids));

	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::ALL_MODELS;


}

TEST_F(RegionManagerTest, BasicTest) {
	Region r = regions[0];
	EXPECT_EQ("r1", r.name);
	EXPECT_EQ(1, r.id);
	EXPECT_EQ(0, r.effStart);
	EXPECT_EQ(25, r.effEnd);
	EXPECT_EQ(10, r.trueStart);
	EXPECT_EQ(15, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ("reg1", r.aliases[0]);
	EXPECT_EQ("1", r.aliases[1]);
	EXPECT_EQ(3, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());

	r = regions[1];
	EXPECT_EQ("r2", r.name);
	EXPECT_EQ(2, r.id);
	EXPECT_EQ(20, r.effStart);
	EXPECT_EQ(50, r.effEnd);
	EXPECT_EQ(30, r.trueStart);
	EXPECT_EQ(40, r.trueEnd);

	r = regions[2];
	EXPECT_EQ("r3", r.name);
	EXPECT_EQ(3, r.id);
	EXPECT_EQ(100, r.effStart);
	EXPECT_EQ(150, r.effEnd);
	EXPECT_EQ(100, r.trueStart);
	EXPECT_EQ(150, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(1, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());
}


TEST_F(RegionManagerTest, TextArchive) {
	BinaryArchive = false;
	regions.WriteArchive("region-test.txt", "\t");
	RegionManager rcopy;
	rcopy.LoadArchive("region-test.txt", "\t");
	Region r = rcopy[0];
	EXPECT_EQ("r1", r.name);
	EXPECT_EQ(0, r.effStart);
	EXPECT_EQ(25, r.effEnd);
	EXPECT_EQ(10, r.trueStart);
	EXPECT_EQ(15, r.trueEnd);
	EXPECT_EQ("reg1", r.aliases[0]);
	EXPECT_EQ("1", r.aliases[1]);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(3, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());

	r = rcopy[1];
	EXPECT_EQ("r2", r.name);
	EXPECT_EQ(20, r.effStart);
	EXPECT_EQ(50, r.effEnd);
	EXPECT_EQ(30, r.trueStart);
	EXPECT_EQ(40, r.trueEnd);

	r = rcopy[2];
	EXPECT_EQ("r3", r.name);
	EXPECT_EQ(100, r.effStart);
	EXPECT_EQ(150, r.effEnd);
	EXPECT_EQ(100, r.trueStart);
	EXPECT_EQ(150, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(1, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());

	remove("region-test.txt");
}


TEST_F(RegionManagerTest, BinaryArchive) {
	BinaryArchive = true;
	regions.WriteArchive("region-test.bin", "\t");
	RegionManager rcopy;
	rcopy.LoadArchive("region-test.bin", "\t");
	Region r = rcopy[0];
	EXPECT_EQ("r1", r.name);
	EXPECT_EQ(0, r.effStart);
	EXPECT_EQ(25, r.effEnd);
	EXPECT_EQ(10, r.trueStart);
	EXPECT_EQ(15, r.trueEnd);
	EXPECT_EQ("reg1", r.aliases[0]);
	EXPECT_EQ("1", r.aliases[1]);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(3, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());

	r = rcopy[1];
	EXPECT_EQ("r2", r.name);
	EXPECT_EQ(20, r.effStart);
	EXPECT_EQ(50, r.effEnd);
	EXPECT_EQ(30, r.trueStart);
	EXPECT_EQ(40, r.trueEnd);

	r = rcopy[2];
	EXPECT_EQ("r3", r.name);
	EXPECT_EQ(100, r.effStart);
	EXPECT_EQ(150, r.effEnd);
	EXPECT_EQ(100, r.trueStart);
	EXPECT_EQ(150, r.trueEnd);
	EXPECT_EQ(2, r.aliases.size());
	EXPECT_EQ(1, r.groups[MetaGroup::DiseaseIndependent].size());
	EXPECT_EQ(2, r.groups[MetaGroup::DiseaseDependent].size());

	remove("region-test.bin");
}

#endif //TEST_APP
