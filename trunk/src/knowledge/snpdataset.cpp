/* 
 * File:   snpdataset.cpp
 * Author: torstees
 * 
 * Created on March 3, 2011, 9:17 AM
 */
#include <iostream>
#include <boost/filesystem.hpp>
#include "snpdataset.h"
#include "utility/exception.h"
#include "utility/strings.h"

namespace Knowledge {

bool SnpDataset::BinaryArchive = false;
uint SnpDataset::rsMapToPositionTolerance = 50;
bool SnpDataset::detailedReport = false;

SnpDataset::SnpDataset(const char *variationsFN) : variationsFilename(variationsFN) {
	markers.reserve(25000);				///< We might have some wasted memory, but in those cases, the requirements will be so small anyway
}

void SnpDataset::SetVariationsFilename(const char *variationsFN) {
	variationsFilename = variationsFN;
}

void SnpDataset::GetFragment(const std::set<std::string>& rsids, SnpDataset& other) {
	Utility::SnpArray::iterator itr = markers.begin();
	Utility::SnpArray::iterator end = markers.end();

	std::set<std::string>::const_iterator missing = rsids.end();
	while (itr != end) {
		if (rsids.find(itr->RSID()) != missing)
			other.AddSNP(itr->chrom, itr->pos, itr->RSID().c_str());
		itr++;
	}
}

void SnpDataset::RsToSnpIndexes(const std::set<std::string>& origrsids, Utility::IdCollection& snps) {
	std::set<std::string>::const_iterator itr = origrsids.begin();
	std::set<std::string>::const_iterator end = origrsids.end();
	std::set<std::string> rsids;			///< These should be upper case to match those in the dataset
	while (itr != end) {
		rsids.insert(Utility::UpperCase(itr++->c_str()));
	}
	
	std::set<std::string>::const_iterator missing = rsids.end();

	uint snpCount = markers.size();

	for (uint i=0; i<snpCount; i++)
		if (rsids.find(markers[i].RSID()) != missing)
			snps.insert(i);
}
void SnpDataset::PositionLookup(Utility::IdCollection& snps, Utility::IdCollection& pos) const {
	Utility::IdCollection::iterator itr = snps.begin();
	Utility::IdCollection::iterator end = snps.end();

	while (itr != end) 
		pos.insert(markers[*itr++].pos);
	
}

void SnpDataset::RangeSnpLookup(char chrom, uint begin, uint end, std::set<LiftOver::SNP>& snps) const {
	Utility::IdCollection indexes;
	chromosomes[(int)chrom - 1].GetSNPs(begin, end, indexes);
	Utility::IdCollection::iterator itr = indexes.begin();
	Utility::IdCollection::iterator iend = indexes.end();

	while (itr != iend) {
		snps.insert(markers[*itr++]);
	}
}

void SnpDataset::RangeSnpLookup(char chrom, uint begin, uint end, Utility::IdCollection& snps) const {
	chromosomes[(int)chrom - 1].GetSNPs(begin, end, snps);
}

void SnpDataset::AddSNPs(const Utility::SnpArray& snps) {
	Utility::SnpArray::const_iterator itr = snps.begin();
	Utility::SnpArray::const_iterator end = snps.end();
	
	while (itr != end) 
		AddSNP(*itr++);
}

uint SnpDataset::AddSNP(const Utility::Locus& l) {
	uint idx = markers.size();
	markers.push_back(l);
	chromosomes[(int)l.chrom - 1].AddSNP(l.pos, idx);
	return idx;
}

uint SnpDataset::AddSNP(char chrom, uint pos, const char *rsid, int role) {
	Utility::Locus snp(chrom, pos, rsid, role);
	//std::cerr<<"Chrom"<<(uint)chrom<<" "<<pos<<" rs"<<rsid<<" "<<role<<"\n";
	uint idx = markers.size();
	markers.push_back(snp);
	chromosomes[(int)chrom - 1].AddSNP(pos, idx);
	return idx;
}

uint SnpDataset::LoadMapData(const char *filename, uint uscBuildVersion, bool performAlignment) {
	std::ifstream file(filename);
	if (!file.good())
		throw Utility::Exception::FileNotFound(filename);
	char line[4096];

	while (file.good()) {
		file.getline(line, 4096);
		std::stringstream ss(line);
		Utility::StringArray words = Utility::Split(line);
		if (words.size() > 3)
			AddSNP(Utility::ChromToInt(words[0].c_str()), atoi(words[3].c_str()), words[1].c_str());
	}

	if (performAlignment)
		AlignData();
	return markers.size();
}

uint SnpDataset::LoadData(const Utility::SnpArray& data, uint uscBuildVersion) {
	// We haven't linked up to liftover yet
	assert(uscBuildVersion==37);

	std::set<Utility::Locus> observedSNPs;		///< To avoid duplicates

	Utility::SnpArray::const_iterator itr = data.begin();
	Utility::SnpArray::const_iterator end = data.end();

	while (itr != end) {
		//Convert marker here!
		if (observedSNPs.find(*itr) == observedSNPs.end()) {
			AddSNP(itr->chrom, itr->pos, itr->RSID().c_str(), itr->role);
			observedSNPs.insert(*itr);
		}
		itr++;
	}

	return observedSNPs.size();
}

void SnpDataset::Clear() {
	markers.clear();
	for (uint i=0; i<Utility::Locus::MaxChromosomes; i++) {
		chromosomes[i].Clear();
	}
}

uint SnpDataset::Size() const {
	return markers.size();
}

Utility::Locus& SnpDataset::operator[](uint index) {
	assert(index<markers.size());
	return markers[index];
}

bool SnpDataset::GetRegionCoverage(uint chrom, uint point, Utility::IdCollection& regionIdxs) {
	if (chrom > 0)
		return chromosomes[chrom-1].GetSegmentCoverage(point, regionIdxs);
	return false;
}

uint SnpDataset::FindSnp(char chrom, const char *rsid, uint pos, uint tolerance) {
	Utility::IdCollection snps;
	uint begin = pos - tolerance;
	uint end = pos + tolerance;

	if (begin > pos)
		begin = 0;
	if (end < pos)
		begin = (uint)-2;
	RangeSnpLookup(chrom, begin, end, snps);

	Utility::IdCollection::iterator itr = snps.begin();
	Utility::IdCollection::iterator iend = snps.end();

	while (itr != iend) {
		if (markers[*itr].RSID() == rsid)
			return *itr;
		itr++;
	}
	return (uint)-1;
}

void SnpDataset::AddRegion(char chrom, uint begin, uint end, uint regionIdx) {
	chromosomes[(int)chrom - 1].AddSegment(begin, end, regionIdx);
}

uint SnpDataset::FindSnp(char chrom, const char *rsid, uint pos) {
	return FindSnp(chrom, rsid, pos, rsMapToPositionTolerance);
}

uint SnpDataset::ReconcileLiftover(std::multimap<Utility::Locus, Utility::Locus>& converted, std::ostream& os, int tolerance) {
	
	os<<"Tolerant SNPs:\n"
	  <<"Chrom\tPos\tRS\n";
	Utility::SnpArray	notFound									= markers;
	if (!boost::filesystem::exists(boost::filesystem::path(variationsFilename)))
		throw Utility::Exception::FileNotFound(variationsFilename.c_str());
	
	std::ifstream file(variationsFilename.c_str(), std::ios::binary);
	if (!file.good()) 
		throw Utility::Exception::FileIO(variationsFilename.c_str(), "Unable to open file.");

	uint offset											= 0;
	uint varVersion									= 0;
	file.read((char*)&varVersion, 4);

	uint count											= 0;
	while (file.good()) {
		char label[3];
		file.read(label, 2);

		if (!file.eof()) {
			label[2]='\0';
			int snpCount								= 0,
					maxPosition							= 0;
			file.read((char*)&snpCount, 4);
			file.read((char*)&maxPosition, 4);

			char chrom = (char)Utility::ChromToInt(label) - 1;
			Chromosome &c								= chromosomes[(int)chrom];
			c.offset										= offset;
			offset+= maxPosition;

			if (file.good()) {
				for (int i=0; i<snpCount; i++) {
					uint rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);

					// If this works, we are success
					uint index = FindSnp(chrom + 1, MakeRSID(rs).c_str(), pos, 1);
					
					if (index == (uint)-1) {
						// These are within tolerance, but will be noted
						index = FindSnp(chrom +1, MakeRSID(rs).c_str(), pos, tolerance);
						
						if (index != (uint)-1) {
							os<<label<<"\t"<<pos<<"\t"<<rs<<"\n";
						}
					}
					
					if (index != (uint)-1) {
						notFound.erase(notFound.begin()+index);
						markers[index].role = role;
						markers[index].pos  = pos;
						count++;
					}
				}
			}
			
			Utility::SnpArray::iterator saItr = notFound.begin();
			Utility::SnpArray::iterator saEnd = notFound.end();
			os<<"\nMissing SNPs:\n"
				 <<"Chrom]\tPos\tRS\n";
			//Work through the errors and report them
			while (saItr != saEnd) {
				os<<Utility::ChromFromInt(saItr->chrom)<<"\t"<<saItr->pos<<"\t"<<saItr->RSID()<<"\n";
				saItr++;
			}
			
		}
	}
	if (count == 0) {
		throw Utility::Exception::FileIO(variationsFilename.c_str(), "No SNPs matched the loaded dataset. Is the variations file in the correct format?");
	}
	return count;
}

uint SnpDataset::AlignData() {
	if (!boost::filesystem::exists(boost::filesystem::path(variationsFilename)))
		throw Utility::Exception::FileNotFound(variationsFilename.c_str());
	
	std::ifstream file(variationsFilename.c_str(), std::ios::binary);
	if (!file.good()) 
		throw Utility::Exception::FileIO(variationsFilename.c_str(), "Unable to open file.");

	uint offset											= 0;
	uint varVersion									= 0;
	file.read((char*)&varVersion, 4);

	uint count											= 0;
	while (file.good()) {
		char label[3];
		file.read(label, 2);

		if (!file.eof()) {
			label[2]='\0';
			int snpCount								= 0,
					maxPosition							= 0;
			file.read((char*)&snpCount, 4);
			file.read((char*)&maxPosition, 4);

			char chrom = (char)Utility::ChromToInt(label) - 1;
			Chromosome &c								= chromosomes[(int)chrom];
			c.offset										= offset;
			offset+= maxPosition;

			if (file.good()) {
				for (int i=0; i<snpCount; i++) {
					uint rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);

					uint index = FindSnp(chrom + 1, MakeRSID(rs).c_str(), pos);
					if (index < (uint)-1) {
						markers[index].role = role;
						count++;
					}
				}
			}
		}
	}
	if (count == 0) {
		throw Utility::Exception::FileIO(variationsFilename.c_str(), "No SNPs matched the loaded dataset. Is the variations file in the correct format?");
	}
	return count;
}
void SnpDataset::SaveArchiveBinary(const char *filename) {
	std::ofstream file(filename, std::ios::binary);
	file.write((char*)&fileVersion, 4);

	for (uint i=0; i<Utility::Locus::MaxChromosomes; i++) {
		Chromosome &chrom = chromosomes[i];
		char ch[3];
		strcpy(ch, Utility::ChromFromInt(i).c_str());
		file.write(ch, 2);

		std::vector<uint> snps = chrom.GetSnpIndexes();
		uint snpCount = snps.size();
		uint offset = chrom.offset;
		file.write((char*)&snpCount, 4);
		file.write((char*)&offset, 4);
		
		for (uint n=0; n<snpCount; n++) {
			Utility::Locus s = markers[snps[n]];
			uint rsid = ExtractRsInt(s.RSID());
			file.write((char*)&rsid, 4);
			file.write((char*)&(s.pos), 4);
			file.write((char*)&(s.role), 1);
		}

	}
}

uint SnpDataset::LoadData(const std::set<std::string>& rsids, std::set<std::string>& snpsFound) {
	return LoadData(variationsFilename.c_str(), rsids, snpsFound);
}

uint SnpDataset::LoadData(const char *filename, const std::set<std::string>& origrsids, std::set<std::string>& snpsFound) {
	if (!boost::filesystem::exists(boost::filesystem::path(filename)))
		throw Utility::Exception::FileNotFound(filename);
	std::ifstream file(filename, std::ios::binary);
	if (!file.good()) 
		throw Utility::Exception::FileIO(filename, "Unable to open file.");

	
	std::set<std::string>::const_iterator itr = origrsids.begin();
	std::set<std::string>::const_iterator end = origrsids.end();
	std::set<std::string> rsids;			///< These should be upper case to match those in the dataset
	while (itr != end) {
		rsids.insert(Utility::UpperCase(itr++->c_str()));
	}
	
	
	
	uint offset											= 0;
	file.read((char*)&fileVersion, 4);

	uint count											= 0;
	while (file.good()) {
		char label[3];
		file.read(label, 2);

		uint localSnps = 0;
		if (!file.eof()) {
			label[2]='\0';
			int snpCount								= 0,
					maxPosition							= 0;
			file.read((char*)&snpCount, 4);
			file.read((char*)&maxPosition, 4);

			char chrom = (char)Utility::ChromToInt(label);
			Chromosome &c								= chromosomes[(int)chrom - 1];
			c.offset										= offset;
			offset+= maxPosition;

			if (file.good()) {
				for (int i=0; i<snpCount; i++) {
					int rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);
					std::string rsid = Utility::_RSID(rs);
					if (rsids.find(rsid) != rsids.end() || rsids.size() == 0) {
						AddSNP(chrom, pos, rsid.c_str(), role);
						snpsFound.insert(rsid);
						localSnps++;
					}
				}
			}
			std::cerr<<"                             Chr "<<label<<" : "<<localSnps<<" SNPs\n";
			count += localSnps;
		}
	}
	std::cerr<<"                  Total SNPs loaded : "<<count<<"\n";
	if (count == 0) 
		throw Utility::Exception::FileIO(filename, "Improper file type");

	return count;
}

void SnpDataset::RoleDescription(uint id, const char *desc) {
	roleDescription[id] = desc;
	std::set<std::string> roles = Utility::ToSet<std::string>(desc, ",");
	std::set<std::string> common;
	std::set_intersection(roles.begin(), roles.end(), exonTerms.begin(), exonTerms.end(), std::inserter(common, common.begin()));
	if (common.size() > 0)
		Utility::Locus::ExonRoles.insert(id);
	common.clear();
	
	std::set_intersection(roles.begin(), roles.end(), intronTerms.begin(), intronTerms.end(), std::inserter(common, common.begin()));
	if (common.size() > 0)
		Utility::Locus::IntronRoles.insert(id);
	common.clear();
	
	std::set_intersection(roles.begin(), roles.end(), regulatoryTerms.begin(), regulatoryTerms.end(),std::inserter(common, common.begin()));
	if (common.size() > 0)
		Utility::Locus::RegulatoryRoles.insert(id);
	common.clear();
	
}

void SnpDataset::LoadFromString(const char *input) {
	if (strlen(input) > 0) {
		Utility::StringArray data=Utility::Split(input, "\t \n", true);
		char chrom	= (char)atoi(data[1].c_str());
		uint pos		= atoi(data[2].c_str());
		uint role	= atoi(data[3].c_str());
		AddSNP(chrom, pos, data[0].c_str(), role);
	}
}

void SnpDataset::LoadArchive(const char *filename) {
	if (BinaryArchive)
		LoadArchiveBinary(filename);
	else {
		char line[1024];
		std::ifstream file(filename);

		while (file.good() && !file.eof()) {
			file.getline(line, 1024);
			LoadFromString(line);
		}
	}
}

void SnpDataset::WriteMarkerInfo(const char *filename, char sep) {
	std::ofstream file(filename);

	uint count = markers.size();
	for (uint i=0; i<count; i++) {
		Utility::Locus &s = markers[i];
		file<<s.RSID()<<sep<<(int)s.chrom<<sep<<s.pos;
		if (detailedReport)
			file<<sep<<roleDescription[s.role];
		file<<"\n";
	}
}


void SnpDataset::SaveArchive(const char *filename, bool writeRole) {
	char sep = '\t';
	if (BinaryArchive)
		SaveArchiveBinary(filename);
	else {
		std::ofstream file(filename);

		uint count = markers.size();
		for (uint i=0; i<count; i++) {
			Utility::Locus &s = markers[i];
			if (writeRole)
				file<<s.RSID()<<sep<<(int)s.chrom<<sep<<s.pos<<sep<<roleDescription[s.role]<<"\n";
			else
				file<<s.RSID()<<sep<<(int)s.chrom<<sep<<s.pos<<sep<<s.role<<"\n";
		}
	/*
		SnpArray::iterator itr = markers.begin();
		SnpArray::iterator end = markers.end();


		uint i =1;
		while (itr != end) {
			std::cout<<i++; std::cout.flush();
			file<<itr->RSID()<<sep<<(int)itr->chr<<sep<<itr->pos<<sep<<itr->role<<"\n";
			std::cout<<"/"<<markers.size()<<" "<<sep<<itr->RSID()<<sep<<(int)itr->chr<<sep<<itr->pos<<sep<<itr->role<<"\n";
			itr++;
		}
	 */
	}
}

void SnpDataset::LoadArchiveBinary(const char *filename) {
	std::set<std::string> rsids;
	std::set<std::string> snpsFound;
	LoadData(filename, rsids, snpsFound);

}




}


#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

TEST(SnpTest, BasicSnpFn) {
	SnpDataset::SNP snp1;
	SnpDataset::SNP snp2(1,5,"rs123",1);
	SnpDataset::SNP snp3(1,10,"rs125",1);
	SnpDataset::SNP snp4(2,1,"rs234",1);

	ASSERT_TRUE(snp2<snp3);
	ASSERT_FALSE(snp3<snp2);
	ASSERT_FALSE(snp4<snp3);
	ASSERT_FALSE(snp2<snp2);

	EXPECT_EQ((uint)5,snp2.Distance(snp3));
	EXPECT_EQ((uint)-1, snp2.Distance(snp4));
	EXPECT_EQ((char)-1, snp1.chrom);
	EXPECT_EQ((uint)-1, snp1.pos);
	EXPECT_EQ("", snp1.RSID());
	EXPECT_EQ(0,snp1.role);

	EXPECT_EQ(1,snp3.chrom);
	EXPECT_EQ(10,snp3.pos);
	EXPECT_EQ("RS125",snp3.RSID());
	EXPECT_EQ(1,snp3.role);
}

class SnpDatasetTest : public ::testing::Test {
public:
	Knowledge::SnpDataset data;
	
	SnpDatasetTest() : data("test.snps"){}

	virtual void SetUp() {
		data.AddSNP(1,1,"rs1",1);		// 0
		data.AddSNP(1,4,"rs10",1);		// 1
		data.AddSNP(1,5,"rs15",1);		// 2
		data.AddSNP(1,10,"rs16",1);	// 3
		data.AddSNP(1,14,"rs20",2);	// 4
		data.AddSNP(1,15,"rs19",2);	// 5
	}
	
	virtual void TearDown() { }
	
};

TEST_F(SnpDatasetTest, TestContents) {
	// Test that there are the correct number inside
	EXPECT_EQ(6,data.Size());
	
	// Use the [] operator to verify the contents
	EXPECT_EQ("RS1", data[0].RSID());
	EXPECT_EQ(1, data[0].pos);

	EXPECT_EQ("RS10", data[1].RSID());
	EXPECT_EQ(4, data[1].pos);
	EXPECT_EQ("RS20", data[4].RSID());
	EXPECT_EQ(14, data[4].pos);
	EXPECT_EQ("RS19", data[5].RSID());
	EXPECT_EQ(15, data[5].pos);
}

TEST_F(SnpDatasetTest, TestQueries) {
	Utility::IdCollection snps;

	std::set<std::string> rsNumbers;
	rsNumbers.insert("rs1");
	rsNumbers.insert("rs10");
	rsNumbers.insert("rs19");
	
	//Let's throw in a couple that aren't present
	rsNumbers.insert("rs14");
	rsNumbers.insert("rs21");

	data.RsToSnpIndexes(rsNumbers, snps);
	EXPECT_EQ((uint)3,snps.size());
	Utility::IdCollection::iterator itr = snps.begin();
	EXPECT_EQ(0,*itr++);
	EXPECT_EQ(1,*itr++);
	EXPECT_EQ(5,*itr++);

	snps.clear();
	data.RangeSnpLookup(1, 3, 10, snps);
	EXPECT_EQ((uint)3, snps.size());
	itr = snps.begin();


	//Make sure that the SNP set isn't being reset when we call one of the functions
	data.RangeSnpLookup(1, 15, 20, snps);
	EXPECT_EQ((uint)4, snps.size());
	EXPECT_EQ(1, *itr++);
	EXPECT_EQ(2, *itr++);
	EXPECT_EQ(3, *itr++);
	EXPECT_EQ(5, *itr);

}

TEST_F(SnpDatasetTest, TestReset) {
	Utility::IdCollection snps;
	data.Clear();
	EXPECT_EQ((uint)0,data.Size());

	data.AddSNP(1,1,"rs1",1);		// 0
	data.AddSNP(1,4,"rs10",1);		// 1
	data.AddSNP(1,5,"rs15",1);		// 2
	data.AddSNP(1,10,"rs16",1);	// 3
	data.AddSNP(1,14,"rs20",2);	// 4
	data.AddSNP(1,15,"rs19",2);	// 5
	data.AddSNP(2,3,"rs80",1);		// 6
	data.AddSNP(2,4,"rs81",1);		// 7
	data.AddSNP(2,5,"rs13",1);		// 8
	data.AddSNP(2,9,"rs17",1);		// 9
	data.AddSNP(2,14,"rs200",2);	// 10
	data.AddSNP(2,15,"rs191",2);	// 11

	EXPECT_EQ(12, data.Size());
}

TEST_F(SnpDatasetTest, TestChrom2) {
	Utility::IdCollection snps;
	data.Clear();

	data.AddSNP(1,1,"rs1",1);		// 0
	data.AddSNP(1,5,"rs10",1);		// 1
	data.AddSNP(1,9,"rs15",1);		// 2
	data.AddSNP(1,10,"rs16",1);	// 3
	data.AddSNP(1,14,"rs20",2);	// 4
	data.AddSNP(1,15,"rs19",2);	// 5
	data.AddSNP(2,3,"rs80",1);		// 6
	data.AddSNP(2,4,"rs81",1);		// 7
	data.AddSNP(2,5,"rs13",1);		// 8
	data.AddSNP(2,9,"rs17",1);		// 9
	data.AddSNP(2,14,"rs200",2);	// 10
	data.AddSNP(2,15,"rs191",2);	// 11

	data.RangeSnpLookup(2,0,5,snps);
	EXPECT_EQ(3,snps.size());
	Utility::IdCollection::iterator itr = snps.begin();
	EXPECT_EQ(6, *itr++);
	EXPECT_EQ(7, *itr++);
	EXPECT_EQ(8, *itr++);
}

void AddChromosome(std::ofstream& file, char *ch, int count, int max) {
	file.write(ch, 2);
	file.write((char*)&count, 4);
	file.write((char*)&max, 4);
}

void AddSNP(std::ofstream& file, int rs, int pos, char role) {
	file.write((char*)&rs, 4);
	file.write((char*)&pos, 4);
	file.write(&role, 1);
}

TEST_F(SnpDatasetTest, TestLoadBadFilename) {
	// Make sure that the file doesn't happen to exist
	remove("test.snps");

	std::set<std::string> rsids;			// We should make sure some aren't found, and don't cover all that are in the file
	rsids.insert("rs1");
	rsids.insert("rs3");
	rsids.insert("rs4");
	rsids.insert("rs7");
	rsids.insert("rs9");
	rsids.insert("rs10");
	rsids.insert("rs20");

	std::set<std::string> found;
	data.Clear();
	// Missing file
	ASSERT_THROW(data.LoadData(rsids, found), Utility::Exception::FileNotFound);
	std::ofstream file("test.snps", std::ios::binary);
	file<<1;
	file.close();

	// Empty file
	ASSERT_THROW(data.LoadData(rsids, found), Utility::Exception::FileIO);

	file.open("test.snps");
	file<<"ABCDEFGHIJ";

	// Meaningless file
	ASSERT_THROW(data.LoadData(rsids, found), Utility::Exception::FileIO);
	remove("test.snps");

}
TEST_F(SnpDatasetTest, TestBinaryAlign) {
	SnpDataset::BinaryArchive = true;
	data.SaveArchive("snps.bin");
	SnpDataset archive("snps.bin");

	SnpDataset::rsMapToPositionTolerance = 5;

	archive.AddSNP(1,1,"rs1");		// 0
	archive.AddSNP(1,6,"rs10");		// 1
	archive.AddSNP(1,9,"rs15");		// 2
	archive.AddSNP(1,10,"rs16");	// 3
	archive.AddSNP(1,13,"rs20");	// 4
	archive.AddSNP(1,20,"rs19");	// 5  (fails to match due to distance)
	archive.AddSNP(2,3,"rx80");		// 6 (fails due to misspelled RS ID)
	archive.AddSNP(2,5,"rs81");		// 7
	archive.AddSNP(2,6,"rs13");		// 8
	archive.AddSNP(2,9,"17");		// 9
	archive.AddSNP(2,14,"rs 200");	// 10 (Fails due to misspelled)
	archive.AddSNP(2,150,"rs191");	// 11

	uint count = archive.AlignData();
	EXPECT_EQ(6, count);

	// Test that there are the correct number inside
	EXPECT_EQ(12, archive.Size());


	// Use the [] operator to verify the contents
	EXPECT_EQ("RS1", archive[0].RSID());
	EXPECT_EQ(1, archive[0].pos);

	EXPECT_EQ("RS10", archive[1].RSID());
	EXPECT_EQ(6, archive[1].pos);
	EXPECT_EQ(1, archive[1].role);
	EXPECT_EQ("RS20", archive[4].RSID());
	EXPECT_EQ(13, archive[4].pos);
	EXPECT_EQ(2, archive[4].role);
	EXPECT_EQ(20, archive[5].pos);
	EXPECT_EQ(2, archive[5].role);
	EXPECT_EQ("RS81", archive[7].RSID());
	EXPECT_EQ(5, archive[7].pos);
	EXPECT_EQ(0, archive[7].role);

	remove("snps.bin");
}
TEST_F(SnpDatasetTest, TestBinaryArchive) {
	SnpDataset::BinaryArchive = true;
	data.SaveArchive("snps.bin");
	SnpDataset archive("");
	archive.LoadArchive("snps.bin");
	
	// Test that there are the correct number inside
	EXPECT_EQ(6, archive.Size());

	// Use the [] operator to verify the contents
	EXPECT_EQ("RS1", archive[0].RSID());
	EXPECT_EQ(1, archive[0].pos);

	EXPECT_EQ("RS10", archive[1].RSID());
	EXPECT_EQ(4, archive[1].pos);
	EXPECT_EQ("RS20", archive[4].RSID());
	EXPECT_EQ(14, archive[4].pos);
	EXPECT_EQ("RS19", archive[5].RSID());
	EXPECT_EQ(15, archive[5].pos);

	remove("snps.bin");
}

TEST_F(SnpDatasetTest, TestTextArchive) {
	SnpDataset::BinaryArchive = false;
	data.SaveArchive("snps.txt");
	SnpDataset archive("");
	archive.LoadArchive("snps.txt");

	// Test that there are the correct number inside
	EXPECT_EQ(6, archive.Size());

	// Use the [] operator to verify the contents
	EXPECT_EQ("RS1", archive[0].RSID());
	EXPECT_EQ(1, archive[0].pos);

	EXPECT_EQ("RS10", archive[1].RSID());
	EXPECT_EQ(4, archive[1].pos);
	EXPECT_EQ("RS20", archive[4].RSID());
	EXPECT_EQ(14, archive[4].pos);
	EXPECT_EQ("RS19", archive[5].RSID());
	EXPECT_EQ(15, archive[5].pos);

	remove("snps.txt");
}

TEST_F(SnpDatasetTest, TestLoadDataFromArray) {
	data.Clear();

	Utility::IdCollection snps;
	SnpDataset::SnpArray snpdata;
	snpdata.push_back(SnpDataset::SNP(1, 1, "rs1", 1));
	snpdata.push_back(SnpDataset::SNP(1, 4, "rs10", 1));
	snpdata.push_back(SnpDataset::SNP(1, 5, "rs15", 1));
	snpdata.push_back(SnpDataset::SNP(1, 10, "rs16", 1));
	snpdata.push_back(SnpDataset::SNP(1, 14, "rs20", 2));
	snpdata.push_back(SnpDataset::SNP(1, 15, "rs19", 2));

	data.LoadData(snpdata, 37);

	EXPECT_EQ(6,data.Size());

	// Use the [] operator to verify the contents
	EXPECT_EQ("RS1", data[0].RSID());
	EXPECT_EQ(1, data[0].pos);

	EXPECT_EQ("RS10", data[1].RSID());
	EXPECT_EQ(4, data[1].pos);
	EXPECT_EQ("RS20", data[4].RSID());
	EXPECT_EQ(14, data[4].pos);
	EXPECT_EQ("RS19", data[5].RSID());
	EXPECT_EQ(15, data[5].pos);
}

TEST_F(SnpDatasetTest, TestLoad) {
	//We are going to create a small variations file to test with
	std::ofstream file("test.snps", std::ios::binary);
	float var = VARIATION_FORMAT_VERSION;
	char chrom[] = "1 ";
	file.write((char*)&var, 4);
	AddChromosome(file, chrom, 5, 15);
	AddSNP(file, 1, 1, 1);		// 0
	AddSNP(file, 2, 5, 1);
	AddSNP(file, 3, 9, 1);		// 1
	AddSNP(file, 4, 12, 1);		// 2
	AddSNP(file, 5, 15, 2);

	chrom[0] = '2';
	AddChromosome(file, chrom, 5, 40);
	AddSNP(file, 3, 3, 1);		// 3
	AddSNP(file, 7, 7, 2);		// 4
	AddSNP(file, 8, 11, 1);
	AddSNP(file, 19, 20, 1);
	AddSNP(file, 20, 22, 1);	// 5

	file.close();
	std::set<std::string> rsids;			// We should make sure some aren't found, and don't cover all that are in the file
	rsids.insert("rs1");
	rsids.insert("rs3");
	rsids.insert("rs4");
	rsids.insert("rs7");
	rsids.insert("rs9");
	rsids.insert("rs10");
	rsids.insert("rs20");

	std::set<std::string> found;
	data.Clear();
	data.LoadData(rsids, found);
	EXPECT_EQ(5, found.size());
	std::set<std::string>::iterator itr = found.begin();
	EXPECT_EQ("RS1", *itr++);
	EXPECT_EQ("RS20", *itr++);
	EXPECT_EQ("RS3", *itr++);
	EXPECT_EQ("RS4", *itr++);
	EXPECT_EQ("RS7", *itr);

	// We should have 6 snps in the dataset, since we have 2 SNPs with RS ID of 3
	EXPECT_EQ(6, data.Size());
	rsids.clear();

	rsids.insert("rs1");
	rsids.insert("rs3");
	rsids.insert("rs20");
	
	
	Utility::IdCollection snps;
	data.RsToSnpIndexes(rsids, snps);
	EXPECT_EQ(4, snps.size());
	Utility::IdCollection::iterator sitr = snps.begin();
	EXPECT_EQ(0, *sitr++);
	EXPECT_EQ(1, *sitr++);
	EXPECT_EQ(3, *sitr++);
	EXPECT_EQ(5, *sitr);

	snps.clear();
	data.RangeSnpLookup(1, 1, 6, snps);
	EXPECT_EQ(1, snps.size());
	EXPECT_EQ(0, *(snps.begin()));

	snps.clear();
	data.RangeSnpLookup(2, 3,7, snps);
	EXPECT_EQ(2, snps.size());
	sitr = snps.begin();
	EXPECT_EQ(3, *sitr++);
	EXPECT_EQ(4, *sitr);

	//Make sure that the lookup doesn't clear pre-existing stuff
	data.RangeSnpLookup(1,10,15,snps);
	sitr = snps.begin();
	EXPECT_EQ(3, snps.size());
	EXPECT_EQ(2, *sitr++);
	EXPECT_EQ(3, *sitr++);
	EXPECT_EQ(4, *sitr);


	//clean up
	remove("test.snps");
}



#endif //TEST_APP
