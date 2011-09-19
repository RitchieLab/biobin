/* 
 * File:   genegenemodelarchive.cpp
 * Author: torstees
 * 
 * Created on March 11, 2011, 1:29 PM
 */

#include "genegenemodelarchive.h"
#include "utility/strings.h"
#include "regionmanager.h"

namespace Knowledge {
	
uint GeneGeneModelArchive::minImplicationIndex = 1.0;		///< Default to all
uint GeneGeneModelArchive::maxModelCount = 10000000;		///< Default to 10 million

void GeneGeneModelArchive::LoadFromArchive(const char *filename, bool useBinary) {
	Reset();
	if (useBinary) {
		std::ifstream file(filename, std::ios::binary);
		uint count = 0;
		file.read((char*)&count, 4);
		while (count-- > 0)
			LoadFromArchive(file, true);
	} else {
		std::ifstream file(filename);
		while (file.good() && !file.eof()) {
			LoadFromArchive(file, false);
		}
	}
}

void GeneGeneModelArchive::GenerateModels(SnpSnpModel::Collection& snpBasedModels, RegionManager& regions) {
	std::set<GeneGeneModel>::iterator itr = models.begin();
	std::set<GeneGeneModel>::iterator end = models.end();

	while (itr != end && snpBasedModels.size() < maxModelCount) 
		itr++->GenerateModels(snpBasedModels, regions);
	
}

void GeneGeneModelArchive::LoadFromArchive(std::istream& os, bool useBinary) {
	GeneGeneModel model;
	if (model.Load(os, useBinary))
		models.insert(model);
}

void GeneGeneModelArchive::WriteToArchive(const char *filename, bool useBinary) {
	if (useBinary) {
		std::ofstream file(filename, std::ios::binary);
		uint count = models.size();
		file.write((char*)&count, 4);
		WriteToArchive(file, true);
	} else {
		std::ofstream file(filename);
		WriteToArchive(file, false);
	}
}


void GeneGeneModelArchive::WriteToArchive(const char *filename, const RegionManager& regions, std::map<float, uint>& scores, bool useBinary) {
	std::ofstream file(filename);
	if (useBinary) {
		file.close();
		file.open(filename, std::ios::binary);
		uint count = models.size();
		file.write((char*)&count, 4);
	}

	iterator itr = models.begin();
	iterator end = models.end();

	uint estimatedModelCount = 0;

	while (itr != end && (estimatedModelCount == (uint)-1 || estimatedModelCount < maxModelCount)) {
		uint n = itr->EstimateModelCount(regions);
		estimatedModelCount += n;
		scores[itr->implicationIndex]+= n;
		//std::cerr<<itr->ToString(regions, "\t")<<"\n";
		itr++->Write(file, useBinary);
	}
}

void GeneGeneModelArchive::WriteToArchive(std::ostream& os, bool useBinary) {
	iterator itr = models.begin();
	iterator end = models.end();

	while (itr != end) {
		itr++->Write(os, useBinary);
	}
}

uint GeneGeneModelArchive::SummarizeModelCounts(std::map<float, uint>& scores, const RegionManager& regions) {
	iterator itr = models.begin();
	iterator end = models.end();

	uint count = 0;
	while (itr != end) {
		uint n = itr->EstimateModelCount(regions);
		count += n;
		scores[itr->implicationIndex]+= n;
		itr++;
	}
	return count;
}

void GeneGeneModelArchive::AddRegions(Utility::IdCollection& ids, RegionManager& regions) {

	if (regions.DoGenerateModels(ids)) {
		Utility::IdCollection visited;
		Utility::IdCollection right = ids;

		Utility::IdCollection::iterator lItr = ids.begin();
		Utility::IdCollection::iterator lEnd = ids.end();

		while (lItr != lEnd) {
			uint lIndex = *lItr;

			if (lIndex < regions.Size()) {
				right.erase(lIndex);
				Utility::IdCollection::iterator rItr = right.begin();
				Utility::IdCollection::iterator rEnd = right.end();

				while (rItr != rEnd) {
					uint rIndex = *rItr;
					if (rIndex < regions.Size()) {
						//Check with the region to validate the model's validity based on configuration
						if (regions.ValidGeneGene(lIndex, rIndex)) {
							float ii = regions[lIndex].ImplicationIndex(regions[rIndex]);
							models.insert(GeneGeneModel(lIndex, rIndex, ii));
						}
					}
					rItr++;
				}
			}
			lItr++;
		}
	}
}

}




#ifdef TEST_APP

#include <gtest/gtest.h>
#include <soci-sqlite3.h>

using namespace Knowledge;

class GeneGeneArchiveTest : public ::testing::Test {
public:
	Knowledge::GeneGeneModelArchive models;
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


		regions.LoadFromDB(sociDB, 1);
		regions.AssociateSNPs(dataset);


		//Lets add some meta IDs so that the implication index won't be 0.0
		Utility::IdCollection meta = Utility::ToSet<uint>("1,2,3,4,5", ",");
		regions[0].AddMetaIDs(MetaGroup::DiseaseIndependent, meta);
		meta = Utility::ToSet<uint>("6", ",");
		regions[0].AddMetaIDs(MetaGroup::DiseaseDependent, meta);

		meta = Utility::ToSet<uint>("1,2,3,4,5", ",");
		regions[1].AddMetaIDs(MetaGroup::DiseaseIndependent, meta);

		meta = Utility::ToSet<uint>("1,2",",");
		regions[2].AddMetaIDs(MetaGroup::DiseaseIndependent, meta);
		meta = Utility::ToSet<uint>("6");
		regions[2].AddMetaIDs(MetaGroup::DiseaseDependent, meta);

		meta = Utility::ToSet<uint>("1,2,4,5", ",");
		regions[3].AddMetaIDs(MetaGroup::DiseaseIndependent, meta);


		//Default II:							Estimated Models
		// 0x1			5.0						12
		// 0x2			3.0						21
		// 0x3			4.0						18
		// 1x2			2.0						28
		// 1x3			4.0						24
		// 2x3			2.0						42

	}

	virtual void TearDown() {
		remove("fake-bio.db");
	}

};
		// Model			II						Est. Model			DD_ONLY?
		// 0x1			5.0						12                 X
		// 0x3			4.0						18                 X
		// 1x3			4.0						24                 0
		// 0x2			3.0						21                 X
		// 1x2			2.0						28                 X
		// 2x3			2.0						42                 X



TEST_F(GeneGeneArchiveTest, Initialization) {
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::ALL_MODELS;

	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,4", ",");
	models.AddRegions(ids, regions);
	EXPECT_EQ(3,models.Size());

	models.Reset();
	EXPECT_EQ(0, models.Size());
	ids = Utility::ToSet<uint>("1,2,0,3", ",");
	models.AddRegions(ids, regions);


	GeneGeneModel test(0,1,5.0);
	GeneGeneModelArchive::iterator itr = models.Begin();
	GeneGeneModel t(*itr++);
	EXPECT_EQ(test, t);
	itr++;			// 0,3, 4.0
	itr++;			// 1,3, 4.0
	test = GeneGeneModel(0,2,3.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(1,2,2.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(2,3,2.0);
	EXPECT_EQ(test, *itr);



}

TEST_F(GeneGeneArchiveTest, CopyConstructor) {
	Utility::IdCollection ids = Utility::ToSet<uint>("0,1,3", ",");
	models.AddRegions(ids, regions);
	EXPECT_EQ(3,models.Size());
	GeneGeneModelArchive newarchive(models);

	GeneGeneModelArchive::iterator i1 = models.Begin();
	GeneGeneModelArchive::iterator i2 = newarchive.Begin();

	EXPECT_EQ(0, (*i1)[0]);
	EXPECT_EQ(1, (*i1)[1]);
	EXPECT_EQ(0, (*i2)[0]);
	EXPECT_EQ(1, (*i2)[1]);
	i2++;
	EXPECT_EQ(0, (*i2)[0]);
	EXPECT_EQ(3, (*i2)[1]);
	i2++;
	EXPECT_EQ(1, (*i2)[0]);
	EXPECT_EQ(3, (*i2)[1]);
	//EXPECT_EQ(*i1++, *i2++);
	//EXPECT_EQ(*i1, *i2);
	
}
		//Default II:							Estimated Models
		// 0x1			5.0						12
		// 0x3			4.0						18
		// 1x2			4.0						28
		// 1x3			4.0						24
		// 0x2			3.0						21
		// 2x3			2.0						42
		// ----------------------------------------
		//                                 145

TEST_F(GeneGeneArchiveTest, ModelSummary) {
	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,3", ",");
	models.AddRegions(ids, regions);
	EXPECT_EQ(6,models.Size());
	// 5.0 -> 12
	// 4.0 -> 18+28+24
	// 3.0 -> 21
	// 2.0 -> 42
	// 1.0 -> 0
	std::map<float, uint> scores;
	uint total = models.SummarizeModelCounts(scores, regions);
	EXPECT_EQ(145, total);
	EXPECT_EQ(12, scores[5.0]);
	EXPECT_EQ(42, scores[4.0]);
	EXPECT_EQ(21, scores[3.0]);
	EXPECT_EQ(70, scores[2.0]);
	EXPECT_EQ(0, scores[1.0]);
}

TEST_F(GeneGeneArchiveTest, GroupOnly) {
	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,3", ",");
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::GROUP_LEVEL;
	models.AddRegions(ids, regions);
	EXPECT_EQ(6,models.Size());
	// 5.0 -> 12
	// 4.0 -> 18+28+24
	// 3.0 -> 21
	// 2.0 -> 42
	// 1.0 -> 0
	std::map<float, uint> scores;
	uint total = models.SummarizeModelCounts(scores, regions);
	EXPECT_EQ(145, total);
	EXPECT_EQ(12, scores[5.0]);
	EXPECT_EQ(42, scores[4.0]);
	EXPECT_EQ(21, scores[3.0]);
	EXPECT_EQ(70, scores[2.0]);
	EXPECT_EQ(0, scores[1.0]);
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::ALL_MODELS;

}
		//Default II:							Estimated          Models
		// 0x1			5.0						12                 X
		// 0x3			4.0						18                 X
		// 1x3			4.0						24                 0
		// 0x2			3.0						21                 X
		// 1x2			2.0						28                 X
		// 2x3			2.0						42                 X
		// --------------------------------------------------------
		//                                 121                 5
TEST_F(GeneGeneArchiveTest, DDOnly) {
	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,3", ",");
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::DD_ONLY;
	models.AddRegions(ids, regions);
	EXPECT_EQ(5,models.Size());
	// 5.0 -> 12
	// 4.0 -> 18
	// 3.0 -> 21
	// 2.0 -> 28 + 42
	// 1.0 -> 0
	std::map<float, uint> scores;
	uint total = models.SummarizeModelCounts(scores, regions);
	EXPECT_EQ(121, total);
	EXPECT_EQ(12, scores[5.0]);
	EXPECT_EQ(18, scores[4.0]);
	EXPECT_EQ(21, scores[3.0]);
	EXPECT_EQ(70, scores[2.0]);
	EXPECT_EQ(0, scores[1.0]);
	RegionManager::modelGenerationType = Knowledge::ModelGenerationMode::ALL_MODELS;
}

TEST_F(GeneGeneArchiveTest, ArchiveText) {
	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,3", ",");
	models.AddRegions(ids, regions);
	models.WriteToArchive("genegenearchive-test.txt", false);
	GeneGeneModelArchive arch;
	arch.LoadFromArchive("genegenearchive-test.txt", false);
	EXPECT_EQ(6, models.Size());
	EXPECT_EQ(6, arch.Size());
	std::map<float, uint> scores;
	uint total = arch.SummarizeModelCounts(scores, regions);
	EXPECT_EQ(145, total);
	EXPECT_EQ(12, scores[5.0]);
	EXPECT_EQ(42, scores[4.0]);
	EXPECT_EQ(21, scores[3.0]);
	EXPECT_EQ(70, scores[2.0]);
	EXPECT_EQ(0, scores[1.0]);
	GeneGeneModel test(0,1,5.0);
	GeneGeneModelArchive::iterator itr = arch.Begin();
	EXPECT_EQ(test, *itr++);
	itr++;			// 0,3
	itr++;			// 1,3
	test = GeneGeneModel(0,2,3.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(1,2,2.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(2,3,2.0);
	EXPECT_EQ(test, *itr);
	remove("genegenearchive-test.txt");
}
		// 0x1			5.0						12                 X
		// 0x3			4.0						18                 X
		// 1x3			4.0						24                 0
		// 0x2			3.0						21                 X
		// 1x2			2.0						28                 X
		// 2x3			2.0						42                 X


TEST_F(GeneGeneArchiveTest, ArchiveBinary) {
	Utility::IdCollection ids = Utility::ToSet<uint>("1,2,0,3", ",");
	models.AddRegions(ids, regions);
	models.WriteToArchive("genegenearchive-test.bin", true);

	GeneGeneModelArchive arch;
	arch.LoadFromArchive("genegenearchive-test.bin", true);
	EXPECT_EQ(6, arch.Size());
	std::map<float, uint> scores;
	uint total = arch.SummarizeModelCounts(scores, regions);
	EXPECT_EQ(145, total);
	EXPECT_EQ(12, scores[5.0]);
	EXPECT_EQ(42, scores[4.0]);
	EXPECT_EQ(21, scores[3.0]);
	EXPECT_EQ(70, scores[2.0]);
	EXPECT_EQ(0, scores[1.0]);
	GeneGeneModel test(0,1,5.0);
	GeneGeneModelArchive::iterator itr = arch.Begin();
	EXPECT_EQ(test, *itr++);
	itr++;			// 0,3
	itr++;			// 1,3
	test = GeneGeneModel(0,2,3.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(1,2,2.0);
	EXPECT_EQ(test, *itr++);
	test = GeneGeneModel(2,3,2.0);
	EXPECT_EQ(test, *itr);

	remove("genegenearchive-test.bin");
}




#endif // TEST_APP
