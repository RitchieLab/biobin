#include "strings.h"


#ifdef TEST_APP
#include <gtest/gtest.h>

using namespace Utility;

TEST(StringsTest, Split) {
	StringArray strings = Split("Region 1,b9028,asdf,,b9028,", ",");
	EXPECT_EQ(4, strings.size());
	EXPECT_EQ("Region 1", strings[0]);
	EXPECT_EQ("b9028", strings[1]);
	EXPECT_EQ("asdf", strings[2]);
	EXPECT_EQ("b9028", strings[3]);

	strings = Split("Region 1,b9028,asdf,,b9028,", ",", true);
	EXPECT_EQ(6, strings.size());

	strings = Split("Region 1 b9028\tasdf\npnng pnng");
	EXPECT_EQ(6, strings.size());
	EXPECT_EQ("Region", strings[0]);
	EXPECT_EQ("1", strings[1]);
	EXPECT_EQ("b9028", strings[2]);
	EXPECT_EQ("asdf", strings[3]);
	EXPECT_EQ("pnng", strings[4]);
	EXPECT_EQ("pnng", strings[5]);
}

TEST(StringTest, SplitTemplated) {
	std::vector<uint> numbers = Split<uint>("1,2,3,4,5,6,7,8,9", ",");
	EXPECT_EQ(9, numbers.size());
	EXPECT_EQ(1, numbers[0]);
	EXPECT_EQ(2, numbers[1]);
	EXPECT_EQ(3, numbers[2]);
}



TEST(StringTest, ToSet) {
	std::set<uint> numbers = ToSet<uint>("2,1,5,4,3,9,7,8,6", ",");
	EXPECT_EQ(9, numbers.size());
	std::set<uint>::iterator itr = numbers.begin();
	EXPECT_EQ(1, *itr++);
	EXPECT_EQ(2, *itr++);
	EXPECT_EQ(3, *itr++);
	EXPECT_EQ(4, *itr++);
	EXPECT_EQ(5, *itr++);
	EXPECT_EQ(6, *itr++);
	EXPECT_EQ(7, *itr++);
	EXPECT_EQ(8, *itr++);
	EXPECT_EQ(9, *itr++);
}

TEST(StringTest, Join2) {
	StringArray strings = Split("Region 1,b9028,asdf,,b9028,", ",", true);
	std::string s = Join<StringArray>(strings, ",");
	EXPECT_EQ("Region 1,b9028,asdf,,b9028,", s);
	std::vector<uint> numbers;

	for (uint i=0; i<10; i++)
		numbers.push_back(i);

	s = Join<std::vector<uint> >(numbers, "-");
	EXPECT_EQ("0-1-2-3-4-5-6-7-8-9", s);
}

TEST(StringTest, Join) {
	StringArray strings = Split("1,2,3,4,5,6", ",");
	std::string result = Join(strings, ",");
	EXPECT_EQ("1,2,3,4,5,6", result);
	std::string result2 = Join(strings, " ");
	EXPECT_EQ("1 2 3 4 5 6", result2);
	std::string result3 = Join(strings, ".", 2, -1);
	EXPECT_EQ("3.4.5.6", result3);
}


TEST(StringTest, TestSplitString) {
	StringArray strings = StringSplit("This is where the GROUP lies. But, GROUP do we break in more than one GROUP place?", "GROUP");
	EXPECT_EQ(4, strings.size());
	EXPECT_EQ("This is where the ", strings[0]);
	EXPECT_EQ(" lies. But, ", strings[1]);
	EXPECT_EQ(" do we break in more than one ", strings[2]);
	EXPECT_EQ(" place?", strings[3]);
}


TEST(StringTest, LoadContentsTest) {
	StringArray strings = Split("String 1,Another String,And one more! But this one is longer.", ",");
	std::ofstream file("load-strings.test");
	StringArray::iterator itr = strings.begin();
	StringArray::iterator end = strings.end();

	while (itr != end) {
		file<<*itr++<<"\n";
	}
	file.close();

	std::string filecontents = LoadContents("load-strings.test");
	EXPECT_EQ("String 1\nAnother String\nAnd one more! But this one is longer.\n", filecontents);
	strings = Split(filecontents.c_str(), "\n");



	remove("load-strings.test");
}

#endif //TEST_APP
