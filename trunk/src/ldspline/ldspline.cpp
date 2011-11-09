/* 
 * File:   ldspline.cpp
 * Author: torstees
 * 
 * Created on August 30, 2010, 12:16 PM
 *
 * Most boundary details associated with loci represent bp locations (probably all of them)
 * Many of the functions here are just redirections to the appropriate call on a given
 * chromosome, so don't get too excited about any interesting algorithmic details here.
 */

#include "ldspline.h"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace Spline {

LdSpline::LdSpline()  {}
LdSpline::LdSpline(const LdSpline& orig) {}
LdSpline::~LdSpline() {}

void LdSpline::LoadHeadersFromHapmap(const char *chrom, const char* filename) {
	LocusLookup loci(&file, chrom);
	filenames[chrom] = filename;
	loci.LoadHeadersFromHapmap(filename);
	this->loci[chrom] = loci;
}

std::pair<int, int> LdSpline::GetBoundariesDP(const char *chrom, int lowerBound, int upperbound, float value) {
	return loci[chrom].GetRangeBoundariesDP(lowerBound, upperbound, value);
}

std::pair<int, int> LdSpline::GetBoundariesRS(const char *chrom, int lowerBound, int upperbound, float value) {
	return loci[chrom].GetRangeBoundariesRS(lowerBound, upperbound, value);
}

void LdSpline::ReportBounds(const char *chrom, int position, float value, const char *type, std::ostream& os) {
	if (loci.find(chrom) != loci.end()) {
		std::pair<int, int> bounds(position, position);

		if (strcmp(type, "dp") == 0 || strcmp(type, "DP") == 0)
			bounds = loci[chrom].GetLdSplineBoundsDP(position, value);
		else if (strcmp(type, "rs") == 0 || strcmp(type, "RS") == 0)
			bounds = loci[chrom].GetLdSplineBoundsRS(position, value);
		else {
			std::cerr<<"Unknown LD type: "<<type<<"\n";
			abort();
		}

		os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\trs"
			 <<loci[chrom].GetPosToRS(bounds.first)<<"\t"<<bounds.first<<"\trs"
			 <<loci[chrom].GetPosToRS(bounds.second)<<"\t"<<bounds.second<<"\t"<<value<<"\n";
	}
	else {
		os<<chrom<<"\t ?? \t"<<position<<"\t ??\t"<<position<<"\t ??\t"<<position<<"\n";
	}
}
void LdSpline::RunReport(const char *chrom, int position, float value, const char *type, std::ostream& os) {
	if (loci.find(chrom) != loci.end()) {
		std::map<int, float> spline;

		if (strcmp(type, "dp") == 0 || strcmp(type, "DP") == 0)
			spline = loci[chrom].GetLdSplineDP(position, value);
		else if (strcmp(type, "rs") == 0 || strcmp(type, "RS") == 0)
			spline = loci[chrom].GetLdSplineRS(position, value);
		else {
			std::cerr<<"Unknown LD type: "<<type<<"\n";
			abort();
		}
		std::map<int, float>::iterator sitr = spline.begin();
		std::map<int, float>::iterator send = spline.end();

		while (sitr != send) {
			std::string direction  = "down";
			if (position < sitr->first)
				direction = "up";
			os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\trs"<<loci[chrom].GetPosToRS(sitr->first)<<"\t"<<sitr->first<<"\t"<< std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setprecision(4)<<sitr->second<<"\t"<<direction<<"\n";
			sitr++;
		}

		if (spline.size() == 0)
			if (loci[chrom].GetPosToRS(position) > 0)
				os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\t \t \n";
			else
				os<<chrom<<"\t  ??  \t"<<position<<"\t \t \n";
	}
	else {
		os<<chrom<<"\t ?? \t"<<position<<"\t"<<"??"<<"\t"<<" "<<"\n";
	}
}



std::vector<SnpSpline> LdSpline::GetLocusRange(const char *chrom, int start, int stop) {
	std::vector<SnpSpline> splines;
	if (loci.find(chrom) != loci.end()) {
		splines = loci[chrom].GetLocusRange(start, stop);
		std::vector<SnpSpline>::iterator itr = splines.begin();
		std::vector<SnpSpline>::iterator end = splines.end();


	}
	return splines;
}

void LdSpline::Summarize(std::ostream& os, const char *chrom) {
	std::map<std::string, LocusLookup>::iterator litr = loci.begin();
	std::map<std::string, LocusLookup>::iterator lend = loci.end();
	while (litr != lend) {
		if (litr->first == chrom || litr->first == "ALL")
			litr->second.Summarize(os);
		litr++;
	}
}

void LdSpline::RunReport(std::ostream& os) {
	std::map<std::string, LocusLookup>::iterator litr = loci.begin();
	std::map<std::string, LocusLookup>::iterator lend = loci.end();

	os<<"Report\n";
	while (litr != lend) {
		LocusLookup::PosToRS rsLookup = litr->second.GetPosToRS();
		LocusLookup::PosToRS::iterator itr = rsLookup.begin();
		LocusLookup::PosToRS::iterator end = rsLookup.end();

		while (itr != end) {
			std::map<int, float> spline =  litr->second.GetLdSplineDP(itr->first, 0.90);
			std::stringstream st100;
			st100<<" DP 0.90: ";
			std::map<int, float>::iterator sitr = spline.begin();
			std::map<int, float>::iterator send = spline.end();
			while (sitr != send) {
				st100<<"\t"<<sitr->first<<" ("<<sitr->second<<")";
				sitr++;
			}

			spline =  litr->second.GetLdSplineDP(itr->first, 0.01);
			std::stringstream st75;
			st75<<" DP 0.01: ";
			sitr = spline.begin();
			send = spline.end();
			while (sitr != send) {
				st75<<"\t"<<sitr->first<<" ("<<sitr->second<<")";
				sitr++;
			}

			os<<itr->first<<"\t"<<itr->second<<"\n\t"<<st75.str()<<"\n\t"<<st100.str()<<"\n";
			itr++;
		}
		litr++;
	}
}

void LdSpline::SaveToCopyBinary(const char *newFilename) {
	file.flush();
	std::fstream newFile(newFilename,  std::ios::out|std::ios::in|std::ios::binary|std::ios::trunc);

	if (newFile.fail()) {
		std::cerr<<"Unable to open file, "<<newFilename<<"\n";
		abort();
	}

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();

	int count = loci.size();

	newFile.write((char*)&count, 4);
	int64_t offset = 4;

	while (itr != end)
		offset += (8 + itr++->second.LocusCount() * 16);

	// At this point, we should have the right offset for the first chromosome's splines
	itr = loci.begin();
	while (itr != end) {
		offset = itr->second.DumpBinaryHeader(&newFile, offset);
		itr++;
	}

	itr = loci.begin();
	while (itr != end) {
		std::string filename = filenames[itr->first];
		itr->second.DumpBinary(&newFile);
		itr++;
	}

	newFile.flush();
}

void LdSpline::SaveToBinary(const char *filename) {
	std::cerr<<"Opening file: "<<filename<<"\n";
	file.close();
	file.open(filename, std::ios::out|std::ios::in|std::ios::binary|std::ios::trunc);

	if (file.fail()) {
		std::cerr<<"Unable to open file, "<<filename<<"\n";
		abort();
	}

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();

	int count = loci.size();
	file.write((char*)&count, 4);

	int64_t offset = 4;

	while (itr != end)
		offset += (8 + itr++->second.LocusCount() * 16);

	// At this point, we should have the right offset for the first chromosome's splines
	itr = loci.begin();
	while (itr != end) {
		offset = itr->second.DumpBinaryHeader(offset);
		itr++;
	}

	itr = loci.begin();
	while (itr != end) {
		std::string filename = filenames[itr->first];
		itr->second.LoadLdFromHapmap(filename.c_str());
		itr->second.DumpBinary();
		itr->second.Release();
		itr++;
	}

	file.flush();

}

void LdSpline::LoadFromBinary(const char *filename) {
	file.open(filename, std::ios::in|std::ios::out|std::ios::binary);
	int count = 0;
	file.read((char*)&count, 4);

	for (int i=0; i<count; i++) {
		LocusLookup chromosome(&file);
		chromosome.LoadBinaryHeader();
		loci[chromosome.Chromosome().c_str()] = chromosome;
	}

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();
	while (itr != end) 
		itr++->second.LoadBinary();
}

/**
 * Basically build the chromosome lookup structure and the barebones
 * stuff required to load a given chromosome when it's time. If they
 * are eager, they can pass true to loadFullheaders, but that's just silly.
 */
void LdSpline::OpenBinary(const char *filename, bool loadFullHeaders) {
	file.open(filename, std::ios::in|std::ios::out|std::ios::binary);
	int count = 0;
	file.read((char*)&count, 4);

	for (int i=0; i<count; i++) {
		LocusLookup chromosome(&file);
		if (loadFullHeaders) 
			chromosome.LoadBinaryHeader();
		else
			chromosome.SkimBinaryHeader();

		loci[chromosome.Chromosome().c_str()] = chromosome;
	}
}

/**
 * A, data abstraction at it's finest:) or something like that
 */
std::map<std::string, LocusLookup> LdSpline::GetChromosomes() {
	return loci;
}


/**
 * This is pretty much a redirection for the lookup's liftover method.
 * The file format is just:
 *
 * chr bp-loc bp-loc+1
 *
 * Lift over seems to require ranges and the chromosomes must start with the
 * letters, chr
 */
void LdSpline::ExportForLiftOver(const char *bimOrig) {
	std::ofstream file(bimOrig);

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();
	while (itr != end) 
		itr++->second.WriteMapForLiftOver(file);
}


/**
 * This is tricky only because there are two files, the primary (loFilename) and
 * loUnmapped, which contains stuff that didn't make it into the new build.
 *
 * Basically, we will load the contents of loUnmapped into a lookup which will
 * be consulted for each locus prior to pulling from the loFilename file. If it's not
 * present, then we are safe to read from the file, if it is, then we know that this
 * didn't make the cut and we move on. When we do encounter something that is missing,
 * we will set the location to -1. It's much too difficult to actually drop nodes
 * out of the spline file, so forever they will be.
 */
void LdSpline::ImportLiftOver(const char *loFilename, const char *loUnmapped) {
	std::map<std::string, std::set<int> > unmappedPositions;		///< Chr -> position [position...]

	//Create a bunch of sets...multimap makes it annoying to do searches for the resulting contents
	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();


	std::ifstream infile(loUnmapped);
	char line[4096];

	//Set up the unampped lookup
	while (infile.good() && !infile.eof()) {
		infile.getline(line, 4096);
		if (line[0]!= '#') {
			std::stringstream ss(line);
			std::string chr;
			int pos = 0;

			ss>>chr>>pos;

			if (pos > 0) {
				std::cerr<<pos<<" - unmapped\n";
				chr.erase(0,3);
				unmappedPositions[chr].insert(pos);
			}
		}
	}

	infile.close();
	infile.open(loFilename);

	itr = loci.begin();
	while (itr != end) {
		itr++->second.LoadMapFromLiftOver(infile, unmappedPositions[itr->first]);
	}
	
}

/**
 * This is just the master function for loading the data according to the import
 * file contents.
 *
 * File format:
 *
 *		chromosome hapmap-filename
 *
 * These should be whitespace separated and appear one per line
 */
std::string LdSpline::ConvertHaploLD(const char *filename) {
	std::ifstream haplo(filename);

	while (haplo.good() && !haplo.eof()) {
		std::string chrom = "", filename = "";
		haplo>>chrom>>filename;

		if (chrom != "")
			LoadHeadersFromHapmap(chrom.c_str(), filename.c_str());
	}


	std::string newFilename = std::string(Utility::ExtractBaseFilename(filename)) + ".ldspline";
	SaveToBinary(newFilename.c_str());

	return newFilename;
}


}
