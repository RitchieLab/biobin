#include "locuslookup.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "utility/strings.h"
#include <cstring>

namespace Spline {

LocusLookup::LocusLookup(std::fstream* file, const char *chr) : chromosome(chr), file(file) { }

LocusLookup::~LocusLookup() { }

void LocusLookup::AddLocus(int rs, int pos, uint64_t offset) {
	std::map<int, int>::iterator itr = posToIdx.find(pos);
	if (itr == posToIdx.end()) {
		int idx = loci.size();
		loci.push_back(SnpSpline(idx, rs, pos, offset));
		posToIdx[pos] = idx;
	}
}

/**
 * We want to skip over any SNP that is found in unmapped, since that didn't make
 * the cut, so to speak, and won't appear in the file we are reading. This is necessary,
 * since we have no way other than order of appearance to relate new positions with
 * old.
 *
 * To denote any location that was unmapped, we will designate it with a -1 for position
 */
void LocusLookup::LoadMapFromLiftOver(std::istream& infile, std::set<int>& unmapped) {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	std::set<int>::iterator umItr = unmapped.begin();
	std::set<int>::iterator umEnd = unmapped.end();



	int idx = 0;
	while (itr != end) {
		std::string chrom, junk;
		int pos = -1;

		idx++;
		if (unmapped.find(itr->pos) == unmapped.end()) {
			char line[4096];
			infile.getline(line, 4096);
			std::stringstream ss(line);

			ss>>chrom>>pos>>junk;
			chrom.erase(0,3);					// They require you to prefix chromosomes as chr1...which we don't want to use
			assert(chromosome==chrom);
		}
		itr++->pos = pos;
	}
}

/**
 * Liftover seems to expect the two boundaries of the bim file to differ, so....
 * let's add one to it and ignore that column when it comes in!
 * @param os
 */
void LocusLookup::WriteMapForLiftOver(std::ostream& os) {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) {
		os<<"chr"<<chromosome<<"\t"<<itr->pos<<"\t"<<itr->pos+1<<"\n";
		itr++;
	}
}

int LocusLookup::GetSnpCount(){
	return loci.size();
}
int LocusLookup::GetPosToRS(int pos) {
	int idx = GetIndex(pos);
	if (idx > 0) {
		return loci[idx].rs;
	}
	return -1;
}



int LocusLookup::GetIndex(int pos) {
	if (posToIdx.find(pos) != posToIdx.end())
		return posToIdx[pos];
	std::map<int, int>::iterator itr = posToIdx.begin();
	std::map<int, int>::iterator end = posToIdx.end();

	return -1;
}

int LocusLookup::GetRS(int idx) {
	return loci[idx].rs;
}

int LocusLookup::GetPos(int idx) {
	return loci[idx].pos;
}

void LocusLookup::Summarize(std::ostream& os) {
	LoadLocusDetails();
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) 
		itr++->Summarize(os, file, chromosome.c_str());
	
}

void LocusLookup::SkimBinaryHeader() {
	char chrom[] = "       ";
	locusCount = 0;
	file->read(chrom, 4);
	chromosome = LDUtility::StripTrailingWhitespace(chrom);
	file->read((char*)&locusCount, 4);
	headerLoadPosition = file->tellg();
	
	file->ignore(16*locusCount);

}

void LocusLookup::LoadLocusDetails() {
	if (loci.size() == 0) {
		file->seekg(headerLoadPosition);

		for (int i=0; i<locusCount; i++) {
			int rs=0, pos=0;
			uint64_t offset=0;
			file->read((char*)&rs, 4);
			file->read((char*)&pos, 4);
			file->read((char*)&offset, 8);
			AddLocus(rs, pos, offset);
		}
	}
}

void LocusLookup::LoadBinaryHeader() {

	char chrom[] = "       ";
	locusCount = 0;
	file->read(chrom, 4);
	chromosome = LDUtility::StripTrailingWhitespace(chrom);
	file->read((char*)&locusCount, 4);
	headerLoadPosition = file->tellg();
	LoadLocusDetails();

	std::fstream::off_type headerPos = file->tellg();
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) {
		itr++->Preload(file);
	}

	file->seekg(headerPos, std::ios::beg);

}

int64_t LocusLookup::DumpBinaryHeader(std::fstream *file, int64_t offset) {
	int locusCount = loci.size();
	char chrom[] = "      ";
	strncpy(chrom, chromosome.c_str(), 4);
	file->write(chrom, 4);
	file->write((char*)&locusCount, 4);

	std::cerr<<chrom<<"<<\t"<<locusCount<<"\t"<<loci.size()<<"\t"<<file->tellp()<<"\n";

	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();


	while (itr != end) {
		int rs = itr->rs;
		int pos = itr->pos;

		assert((int64_t)itr->offset == offset);
		itr->offset = offset;
		file->write((char*)&rs, 4);
		file->write((char*)&pos, 4);
		file->write((char*)&offset, 8);
		offset += (8 + itr->GetSplineCount() * 8);

		itr++;
	}
	return offset;
}
/**
 * Dumps item to a binary file.
 * Binary format:
 *		* 4 bytes *  Chromosome string (right now, only 1-22 and X & Y)
 *		* 4 bytes *  Number of Loci (total number of loci in file)
 *		* For each locus *
 *			* 4 bytes * rsID
 *			* 4 bytes * pos
 *       * 8 bytes * offset (this is where to jump to to start reading data for that locus
 * @param filename
 */
int64_t LocusLookup::DumpBinaryHeader(int64_t offset) {

	int locusCount = loci.size();
	char chrom[] = "      ";
	strncpy(chrom, chromosome.c_str(), 4);
	file->write(chrom, 4);
	file->write((char*)&locusCount, 4);

	std::cerr<<chrom<<"<<\t"<<locusCount<<"\t"<<file->tellp()<<"\n";

	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	
	while (itr != end) {
		int rs = itr->rs;
		int pos = itr->pos;
		itr->offset = offset;
		file->write((char*)&rs, 4);
		file->write((char*)&pos, 4);
		file->write((char*)&offset, 8);
		offset += (8 + itr->GetSplineCount() * 8);

		itr++;
	}
	return offset;
}

int LocusLookup::LocusCount() {
	return loci.size();
}

void LocusLookup::LoadBinary() {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	int idx = 0;
	while (itr != end) {
		itr->LoadFromBinary(file);
		itr++;
		idx += 1;
	}
}
/**
 * Dumps the ld splines to a binary file
 * Binary Format: (for each SNP, the spline is recorded)
 *		* 4 bytes * count
 *			* 4 bytes * DPrime
 *			* 4 bytes * RSquared
 
 * @param filename
 */
void LocusLookup::DumpBinary() {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr!= end) {

		itr->WriteToBinary(file);
		itr++;
	}

}
void LocusLookup::DumpBinary(std::fstream* file) {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr!= end) {
		bool doRelease = itr->LoadFromBinary(this->file);
		itr->WriteToBinary(file);
		if (doRelease)
			itr->ReleaseStats();
		itr++;
	}

}
std::string LocusLookup::Chromosome() const {
	return chromosome;
}

void LocusLookup::IncrementSplineCounts(int pos1, int pos2) {
	int id1 = GetIndex(pos1);
	int id2 = GetIndex(pos2);

	loci[id1].UpdateUpstreamSplineCount(id2);
	loci[id2].UpdateDownstreamSplineCount(id1);
}

void LocusLookup::AddLdValue(int pos1, int pos2, float dp, float rs) {
	int id1 = GetIndex(pos1);
	int id2 = GetIndex(pos2);

	//For now, we are under the assumption that pos1<pos2
	loci[id1].AddStatsUpstream(id2, dp, rs);
	loci[id2].AddStatsDownstream(id1, dp, rs);

}

LocusLookup::RsToPos LocusLookup::GetRsToPos() {
	RsToPos rsToPos;

	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) {
		rsToPos[itr->rs] = itr->pos;
		itr++;
	}

	return rsToPos;
}


LocusLookup::PosToRS LocusLookup::GetPosToRS() {
	PosToRS posToRS;

	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) {
		//std::cerr<<"  "<<itr->pos<<" -> "<<itr->rs<<"\n";
		posToRS[itr->pos] = itr->rs;
		itr++;
	}

	return posToRS;
}

/**
 * Start and stop are just min/max boundaries for the region we are searcing over. They
 * don't actually have to appear as observed loci in the original hapmap data
 */
std::pair<int, int> LocusLookup::GetRangeBoundariesRS(int start, int stop, float value ) {
	LoadLocusDetails();
	
	std::pair<int, int> bounds(start, stop);



	std::map<int, int>::iterator lb = posToIdx.upper_bound(start-1);
	std::map<int, int>::iterator ub = posToIdx.upper_bound(stop);

	while (lb != ub) {
		SnpSpline &sp = loci[lb->second];
		bool doRelease = sp.LoadFromBinary(file);
		std::pair<int, int> localBounds = sp.GetSplineBoundsRS(value);

		if (loci[localBounds.first].pos > 0 && loci[localBounds.second].pos > 0) {
			if (loci[localBounds.first].pos < bounds.first)
				bounds.first = loci[localBounds.first].pos;
			if (loci[localBounds.second].pos > bounds.second)
				bounds.second = loci[localBounds.second].pos;
		}

		if (doRelease)
			sp.ReleaseStats();
		lb++;
	}
	return bounds;
}

/**
 * Start and stop are just min/max boundaries for the region we are searcing over. They
 * don't actually have to appear as observed loci in the original hapmap data
 */
std::pair<int, int> LocusLookup::GetRangeBoundariesDP(int start, int stop, float value ) {
	LoadLocusDetails();
	std::pair<int, int> bounds(start, stop);

	std::map<int, int>::iterator lb = posToIdx.upper_bound(start-1);
	std::map<int, int>::iterator ub = posToIdx.upper_bound(stop);

	//uint posToIdxSize = posToIdx.size();

	while (lb != ub) {
		SnpSpline &sp = loci[lb->second];
		bool doRelease = sp.LoadFromBinary(file);
		std::pair<int, int> localBounds = sp.GetSplineBoundsDP(value);
		if (loci[localBounds.first].pos > 0 && loci[localBounds.second].pos > 0) {
			if (loci[localBounds.first].pos < bounds.first)
				bounds.first = loci[localBounds.first].pos;
			if (loci[localBounds.second].pos > bounds.second)
				bounds.second = loci[localBounds.second].pos;

		}
		if (doRelease)
			sp.ReleaseStats();
		lb++;
	}


	return bounds;
}

/**
 * Start and stop are just min/max boundaries for the region we are searcing over. They
 * don't actually have to appear as observed loci in the original hapmap data
 */
std::vector<SnpSpline> LocusLookup::GetLocusRange(int start, int stop) {
	LoadLocusDetails();
	std::vector<SnpSpline> stuff;
	std::map<int, int>::iterator lb = posToIdx.upper_bound(start-1);
	std::map<int, int>::iterator ub = posToIdx.upper_bound(stop);

	while (lb != ub) {
		stuff.push_back(loci[lb->second]);
		lb++;
	}
	return stuff;
}

std::pair<int, int> LocusLookup::GetLdSplineBoundsDP(int pos, float value) {
	LoadLocusDetails();

	int idx = GetIndex(pos);

	std::pair<int, int> dp(pos, pos);
	if (idx > 0) {
		SnpSpline &spline = loci[idx];
		bool doRelease = spline.LoadFromBinary(file);

		std::pair<int, int> indexes = spline.GetSplineBoundsDP(value);
		dp.first = loci[indexes.first].pos;
		dp.second = loci[indexes.second].pos;

		if (doRelease)
			spline.ReleaseStats();
	}
	return dp;
}

std::pair<int, int> LocusLookup::GetLdSplineBoundsRS(int pos, float value) {
	LoadLocusDetails();

	int idx = GetIndex(pos);

	std::pair<int, int> rs(pos, pos);
	if (idx > 0) {
		SnpSpline &spline = loci[idx];
		bool doRelease = spline.LoadFromBinary(file);

		std::pair<int, int> indexes = spline.GetSplineBoundsRS(value);
		rs.first = loci[indexes.first].pos;
		rs.second = loci[indexes.second].pos;
		if (doRelease)
			spline.ReleaseStats();
	}
	return rs;
}


LocusLookup::SplineLookup LocusLookup::GetLdSplineDP(int pos, float dp) {
	SplineLookup positions;
	LoadLocusDetails();

	int idx = GetIndex(pos);

	if (idx > 0) {
		SnpSpline &spline = loci[idx];
		bool doRelease = spline.LoadFromBinary(file);
		std::map<int, float> indexes = spline.GetSplineDP(dp);

		std::map<int, float>::iterator itr = indexes.begin();
		std::map<int, float>::iterator end = indexes.end();

		while (itr != end) {
			positions[GetPos(itr->first)] = itr->second;
			itr++;
		}

		if (doRelease)
			spline.ReleaseStats();
	}
	return positions;
}



LocusLookup::SplineLookup LocusLookup::GetLdSplineRS(int pos, float dp) {
	SplineLookup positions;
	LoadLocusDetails();

	int idx = GetIndex(pos);

	if (idx > 0) {
		SnpSpline &spline = loci[idx];
		bool doRelease = spline.LoadFromBinary(file);
		std::map<int, float> indexes = spline.GetSplineRS(dp);

		std::map<int, float>::iterator itr = indexes.begin();
		std::map<int, float>::iterator end = indexes.end();

		while (itr != end) {
			positions[GetPos(itr->first)] = itr->second;
			itr++;
		}

		if (doRelease)
			spline.ReleaseStats();
	}
	return positions;
}

void LocusLookup::LoadHeadersFromHapmap(const char* filename) {
	std::ifstream file(filename);
	int offset = 0;							// chrom count offset
	while (file.good() && !file.eof()) {
		int pos1 = -1, pos2 = -1;
		std::string pop = "", rs1 = "", rs2 = "", junk = "";
		float dp = -1.0, rs = -1.0, lod = -1.0;

		file>>pos1>>pos2>>pop>>rs1>>rs2>>dp>>rs>>lod>>junk;

		if (pos1 > 0) {
			int rsid = LDUtility::ExtractRsNumber(rs1.c_str());
			AddLocus(rsid, pos1, offset);
			rsid = LDUtility::ExtractRsNumber(rs2.c_str());
			AddLocus(rsid, pos2, offset);
			IncrementSplineCounts(pos1, pos2);

		}

	}

	std::cerr<<"Snps Parsed from file, "<<filename<<"("<<this->chromosome<<") : "<<loci.size()<<"\n";
}


void LocusLookup::LoadLdFromHapmap(const char *filename) {
	std::ifstream file(filename);

	while (file.good() && !file.eof()) {

		int pos1 = -1, pos2 = -1;
		std::string pop = "", rs1 = "", rs2 = "", junk = "";
		float dp = -1.0, rs = -1.0, lod = -1.0;

		file>>pos1>>pos2>>pop>>rs1>>rs2>>dp>>rs>>lod>>junk;
		if (pos1 > 0 && pos2 > 0)
			AddLdValue(pos1, pos2, dp,rs);

	}
}

void LocusLookup::Release() {
	std::vector<SnpSpline>::iterator itr = loci.begin();
	std::vector<SnpSpline>::iterator end = loci.end();

	while (itr != end) 
		itr++->ReleaseStats();

}

}
