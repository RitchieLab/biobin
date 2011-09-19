//
// C++ Implementation: biomodel.h
//
// Description: Encapsulates any functionality duplicated in Snp/Snp or Gene/Gene models
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __BIOFILTER_BIOMODEL_H
#define __BIOFILTER_BIOMODEL_H

namespace Knowledge {

class BioModel {
public:
	BioModel();
	BioModel(float implIdx);
	virtual ~BioModel();

	virtual float ImplicationIndex() const;
protected:
	float implicationIndex;
};

inline
BioModel::BioModel() { }

inline
BioModel::BioModel(float implIndex) : implicationIndex(implIndex) { }

inline
BioModel::~BioModel() { }
	
inline
float BioModel::ImplicationIndex() const {
	return implicationIndex;
}

}

#endif

