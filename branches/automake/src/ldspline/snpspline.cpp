#include "snpspline.h"
#include <algorithm>
#include <iostream>
#include <assert.h>

namespace Spline {

SnpSpline::SnpSpline(int idx, int rsid, int pos, uint64_t offset) :pos(pos), rs(rsid), idx(idx), offset(offset), usIdx(0), dsIdx(0), pins(0) {}

SnpSpline::~SnpSpline() {}

void SnpSpline::AddStatsUpstream(int idx, float dp, float rs) {
	if (idx >= 0) {
		int lIdx = idx - this->idx - 1;
		if (lIdx >= (int)upstream.size())
			upstream.resize(lIdx + 1);
		upstream[lIdx] = LdStat(dp, rs);
	}
}

void SnpSpline::AddStatsDownstream(int idx, float dp, float rs) {

	if (idx >= 0) {
		int lIdx = this->idx - idx - 1;
		if (lIdx >= (int)downstream.size())
			downstream.resize(lIdx + 1);
		downstream[lIdx] = LdStat(dp, rs);
	}
}

int SnpSpline::GetSplineCount() {
	if (usIdx + dsIdx == 0)
		 return upstream.size() + downstream.size();
	return usIdx + dsIdx;
}

std::vector<LdStat> SnpSpline::GetUpstream() {
	return upstream;
}

std::vector<LdStat> SnpSpline::GetDownstream() {
	return downstream;
}



void SnpSpline::UpdateUpstreamSplineCount(int idx) {
	int d = idx - this->idx;
	if (d > usIdx)
		usIdx = d;
}

void SnpSpline::UpdateDownstreamSplineCount(int idx) {
	int d = this->idx - idx;
	if (d > dsIdx )
		dsIdx = d;
}


std::map<int, float> SnpSpline::GetSplineRS(float minRS) {
	std::map<int, float> indexes;
	std::vector<LdStat>::iterator itr = downstream.begin();
	std::vector<LdStat>::iterator end = downstream.end();

	int i=-1;
	bool doContinue = true;
	while (itr!=end && doContinue) {
		doContinue = itr->rs == -1 || itr->rs >= minRS;
		if (doContinue) {
			if (itr->rs > 0.0) {
			indexes[this->idx+i] = itr->rs;
			}

		}
		i--;
		itr++;
	}

	doContinue = true;
	itr = upstream.begin();
	end = upstream.end();
	i = 0;
	while (itr!=end && doContinue) {
		i++;
		doContinue = (itr->rs == -1 || itr->rs>=minRS);
		if (doContinue) {
			if (itr->rs > 0){
				indexes[this->idx+i]= itr->rs;
			}
		}
		itr++;
	}
	return indexes;
}

std::pair<int, int> SnpSpline::GetSplineBoundsDP(float minDP) {
	int min = idx, max = idx;
	std::vector<LdStat>::iterator itr = downstream.begin();
	std::vector<LdStat>::iterator end = downstream.end();

	bool doContinue = true;
	while (itr!=end && doContinue) {
		doContinue = itr->dp == -1 || itr->dp >= minDP;
		if (doContinue) {
			if (itr->dp > 0.0)
				min -= 1;
		}
		itr++;
	}

	doContinue = true;
	itr = upstream.begin();
	end = upstream.end();
	while (itr!=end && doContinue) {
		doContinue = (itr->dp == -1 || itr->dp>=minDP);
		if (doContinue) {
			if (itr->dp > 0)
				max+=1;
		}
		itr++;
	}
	return std::pair<int,int>(min,max);
}


std::pair<int, int> SnpSpline::GetSplineBoundsRS(float minDP) {
	int min = idx, max = idx;
	std::vector<LdStat>::iterator itr = downstream.begin();
	std::vector<LdStat>::iterator end = downstream.end();

	bool doContinue = true;
	while (itr!=end && doContinue) {
		doContinue = itr->rs == -1 || itr->rs >= minDP;
		if (doContinue) {
			if (itr->rs > 0.0)
				min -= 1;
		}
		itr++;
	}

	doContinue = true;
	itr = upstream.begin();
	end = upstream.end();
	while (itr!=end && doContinue) {
		doContinue = (itr->rs == -1 || itr->rs>=minDP);
		if (doContinue) {
			if (itr->rs > 0)
				max+=1;
		}
		itr++;
	}
	return std::pair<int,int>(min,max);
}



std::map<int, float> SnpSpline::GetSplineDP(float minDP) {
	std::map<int, float> indexes;
	std::vector<LdStat>::iterator itr = downstream.begin();
	std::vector<LdStat>::iterator end = downstream.end();

	int i=-1;
	bool doContinue = true;
	while (itr!=end && doContinue) {
		doContinue = itr->dp == -1 || itr->dp >= minDP;
		if (doContinue) {
			if (itr->dp > 0.0) {
			indexes[this->idx+i] = itr->dp;
			}

		}
		i--;
		itr++;
	}

	//std::reverse(indexes.begin(), indexes.end());
	doContinue = true;
	itr = upstream.begin();
	end = upstream.end();
	i = 0;
	while (itr!=end && doContinue) {
		i++;
		doContinue = (itr->dp == -1 || itr->dp>=minDP);
		if (doContinue) {
			if (itr->dp > 0){
				indexes[this->idx+i]= itr->dp;
			}
		}
		itr++;
	}
	return indexes;
}

void SnpSpline::Summarize(std::ostream& os, std::fstream* file, const char *chrom) {
	file->seekg(offset, std::ios::beg);
	unsigned     int dscount, uscount;
	file->read((char*)&dscount, 4);
	file->ignore(8*dscount);
	file->read((char*)&uscount, 4);
	os<<idx<<"\t"<<chrom<<"\trs"<<rs<<"\t"<<pos<<"\t"<<offset<<"\t"<<dscount<<"\t"<<uscount<<"\n";
}

void SnpSpline::WriteToBinary(std::fstream *file) {
	unsigned int count = downstream.size();

unsigned int dscount = (unsigned int)downstream.size();
unsigned int uscount = (unsigned int)upstream.size();
assert(offset == file->tellp());
bool willfail = (dscount+uscount) != (unsigned int)GetSplineCount();
	offset		= file->tellp();

	file->write((char*)&count, 4);
	std::vector<LdStat>::iterator sItr = downstream.begin();
	std::vector<LdStat>::iterator sEnd = downstream.end();
	while (sItr != sEnd) {
		file->write((char*)&(sItr->dp), 4);
		file->write((char*)&(sItr->rs), 4);
	
		sItr++;
	}


	count = upstream.size();
	file->write((char*)&count, 4);
	sItr = upstream.begin();
	sEnd = upstream.end();

	while (sItr != sEnd) {
		file->write((char*)&(sItr->dp), 4);
		file->write((char*)&(sItr->rs), 4);
		sItr++;
	}

assert((unsigned int)this->GetSplineCount() == (downstream.size() + upstream.size()));
}

void SnpSpline::ReleaseStats() {
	if (--pins < 1) {
		//Actually reduce the memory consumed by the structures
		std::vector<LdStat> u;
		upstream.swap(u);

		std::vector<LdStat> d;
		downstream.swap(d);
	}
}

void SnpSpline::Preload(std::fstream *file) {
	if (downstream.size() + upstream.size() == 0) {
		file->seekg(offset, std::ios::beg);
		dsIdx = usIdx = -1;
		file->read((char*)&dsIdx, 4);
		file->ignore(8*dsIdx);
		file->read((char*)&usIdx, 4);
		
	}
}

bool SnpSpline::LoadFromBinary(std::fstream *file) {
	bool didLoad = false;
	pins++;
	if (downstream.size() + upstream.size() == 0) {


		downstream.reserve(dsIdx);
		upstream.reserve(usIdx);
		didLoad = true;
		//Downstream part
		int count = -1;
		file->seekg(offset, std::ios::beg);
		file->read((char*)&count, 4);
		for (int i=0; i<count; i++){
			float dp=0, rs=0.0;
			file->read((char*)&dp, 4);
			file->read((char*)&rs, 4);
			AddStatsDownstream(idx - i - 1, dp, rs);
		}
		count = -1;
		file->read((char*)&count, 4);
		for (int i=0; i<count; i++) {
			float dp=0, rs=0.0;
			file->read((char*)&dp, 4);
			file->read((char*)&rs, 4);
			AddStatsUpstream(idx + i + 1, dp, rs);
		}
	}
	return didLoad;
}

}
