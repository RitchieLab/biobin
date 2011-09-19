//
// C++ Implementation: snpsnpmodel.cpp
//
// Description: Defines functionality associated with Snp/Snp models
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "snpsnpmodel.h"
#include <fstream>


#ifdef TEST_APP

#include <gtest/gtest.h>

using namespace Knowledge;

TEST(SnpSnpModelTest, BasicFunctionality) {
	SnpSnpModel model1(1,2,0.5);
	SnpSnpModel model2(3,1,2.1);

	EXPECT_TRUE(model1<model2);
	EXPECT_FALSE(model2<model1);

	SnpSnpModel model1a(model1);
	EXPECT_TRUE(model1a==model1);
	EXPECT_FALSE(model1a==model2);

	model1a = model2;
	EXPECT_FALSE(model1a==model1);

	EXPECT_EQ(1, model1[0]);
	EXPECT_EQ((uint)2, model1.Size());
	EXPECT_EQ(1, model2[0]);

	SnpSnpModel model3(2,2,4.0);
	EXPECT_EQ(1, model3.Size());
}

TEST(SnpSnpModelTest, Collection) {
	SnpSnpModel::Collection models;

	models.insert(SnpSnpModel(1,2,5.0));
	models.insert(SnpSnpModel(3,4,1.0));
	models.insert(SnpSnpModel(2,3,2.0));
	models.insert(SnpSnpModel(1,3,3.0));
	models.insert(SnpSnpModel(1,2,2.0));

	EXPECT_EQ((uint)4, models.size());
	SnpSnpModel::Collection::iterator itr = models.begin();

	EXPECT_EQ(1, (*itr)[0]);
	EXPECT_EQ(2, (*itr)[1]);
	EXPECT_EQ(5.0, itr++->ImplicationIndex());

	EXPECT_EQ(1, (*itr)[0]);
	EXPECT_EQ(3, (*itr)[1]);
	EXPECT_EQ(3.0, itr++->ImplicationIndex());

	EXPECT_EQ(2, (*itr)[0]);
	EXPECT_EQ(3, (*itr)[1]);
	EXPECT_EQ(2.0, itr++->ImplicationIndex());

	EXPECT_EQ(3, (*itr)[0]);
}

TEST(SnpSnpModelTest, BinaryArchive) {
	BinaryArchive = true;
	SnpSnpModel::Collection models;
	models.insert(SnpSnpModel(1,2,5.0));
	models.insert(SnpSnpModel(3,4,1.0));
	models.insert(SnpSnpModel(2,3,2.0));
	models.insert(SnpSnpModel(1,3,3.0));
	models.insert(SnpSnpModel(1,2,2.0));

	SnpSnpModel::Collection::iterator itr = models.begin();
	SnpSnpModel::Collection::iterator end = models.end();

	std::ofstream file("snpmodel.test");
	while (itr != end) {
		itr->Write(file);
		itr++;
	}
	file.close();

	std::ifstream infile("snpmodel.test");
	SnpSnpModel model;

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(1, model[0]);
	EXPECT_EQ(2, model[1]);
	EXPECT_EQ(5.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(1, model[0]);
	EXPECT_EQ(3, model[1]);
	EXPECT_EQ(3.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(2, model[0]);
	EXPECT_EQ(3, model[1]);
	EXPECT_EQ(2.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(3, model[0]);
	EXPECT_EQ(4, model[1]);
	EXPECT_EQ(1.0, model.ImplicationIndex());

	remove("snpmodel.test");
}

TEST(SnpSnpModelTest, TextArchive) {
	BinaryArchive = false;
	SnpSnpModel::Collection models;
	models.insert(SnpSnpModel(1,2,5.0));
	models.insert(SnpSnpModel(3,4,1.0));
	models.insert(SnpSnpModel(2,3,2.0));
	models.insert(SnpSnpModel(1,3,3.0));
	models.insert(SnpSnpModel(1,2,2.0));

	SnpSnpModel::Collection::iterator itr = models.begin();
	SnpSnpModel::Collection::iterator end = models.end();

	std::ofstream file("snpmodel.test");
	while (itr != end) {
		itr->Write(file);
		itr++;
	}
	file.close();

	std::ifstream infile("snpmodel.test");
	SnpSnpModel model;

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(1, model[0]);
	EXPECT_EQ(2, model[1]);
	EXPECT_FLOAT_EQ(5.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(1, model[0]);
	EXPECT_EQ(3, model[1]);
	EXPECT_FLOAT_EQ(3.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(2, model[0]);
	EXPECT_EQ(3, model[1]);
	EXPECT_FLOAT_EQ(2.0, model.ImplicationIndex());

	model.Load(infile);
	EXPECT_EQ(2, model.Size());
	EXPECT_EQ(3, model[0]);
	EXPECT_EQ(4, model[1]);
	EXPECT_FLOAT_EQ(1.0, model.ImplicationIndex());

	remove("snpmodel.test");
}

#endif //TEST_APP
