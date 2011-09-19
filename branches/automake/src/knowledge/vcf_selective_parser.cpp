/* 
 * File:   vcf_selective_parser.cpp
 * Author: torstees
 * 
 * Created on April 4, 2011, 3:43 PM
 */

#include <string>
#include <iostream>
#include <fstream>
#include "vcf_selective_parser.h"
#include "utility/rs.h"
#include "utility/strings.h"
#include "chromosome.h"

#include <vector>
using namespace Knowledge;
using namespace Utility;
typedef std::vector<Chromosome> ChromosomeContainer;

struct RsID {
	RsID(const std::string& chrom, uint pos, const std::string& rsid, char ref, char alt ) : chrom(chrom), rsid(rsid), pos(pos), ref(ref), alt(alt) { }

	std::string chrom;
	std::string rsid;
	uint pos;
	char ref;
	char alt;

	friend std::ostream &operator<<(std::ostream &stream, RsID rs) {
		stream<<rs.rsid;
		return stream;
	}
};

struct Individual {
	Individual(const std::string& name) : id(name) {}
	Individual(const Individual& other) : id(other.id), genotypes(other.genotypes) {}
	std::string id;
	std::vector<char> genotypes;

	void AddGenotype(const std::string& s) {

		char genotype = '.' - 1;
		Utility::StringArray words = Utility::Split(s.c_str(), ":");
		if (words[0] == "0|0")
			genotype = '0';
		else if (words[0] == "0|1" || words[0] == "1|0")
			genotype = '1';
		else if (words[0] == "1|1")
			genotype = '2';
		genotypes.push_back(genotype);
	}

	char operator[](const int idx) {
		return genotypes[idx];
	}

	friend std::ostream &operator<<(std::ostream &stream, Individual ind) {
		stream<<ind.id;
		return stream;
	}
	std::string Report() {
		return id + "\t" + Utility::Join(genotypes, "\t");
	}
};

void ParseRegionDefinition(ChromosomeContainer& chromosomes, const char *filename) {
	std::ifstream file(filename);

	while (file.good()) {
		uint start, stop;
		std::string chrom = "";
		char line[8102];
		file.getline(line, 8102);

		Utility::StringArray strings = Utility::Split(line, "\t ,");
		if (strings.size() > 0) {
			chrom = strings[0];
			start = atoi(strings[1].c_str());
			stop  = atoi(strings[2].c_str());

			//file>>chrom>>start>>stop;
			if (chrom.size() > 0) {
				uint ch = Utility::ChromToInt(chrom.c_str());
				chromosomes[ch].AddSegment(start, stop);
			}
		}
	}


}

void ParseVCF(ChromosomeContainer& chromosomes, const char *filename) {
	std::ifstream file(filename);
	std::vector<Individual> dataset;
	std::vector<RsID> rsids;
	char line[1048576];
	while (file.good()) {
		line[0] = '\0';
		file.getline(line, 1048576);
		Utility::StringArray words = Utility::Split(line);
		uint wordCount = words.size();
		if (wordCount > 0) {
			if (words[0] == "#CHROM") {
				for (uint i=9; i<wordCount; i++) {
					dataset.push_back(Individual(words[i]));
				}
			}
			else if (words[0][0] != '#') {
				uint ch = Utility::ChromToInt(words[0].c_str());
				uint pos = atoi(words[1].c_str());
				char ref = words[3][0];
				char alt = words[4][0];
				if (chromosomes[ch].IsCoveredBySegment(pos)) {
					rsids.push_back(RsID(words[0], pos, words[2], ref, alt));
					for (uint i=9; i<wordCount; i++){
						dataset[i-9].AddGenotype(words[i]);
					}
				}
			}
		}
	}

	std::string snpdata = std::string(filename) + ".genotypes";
	std::ofstream outfile(snpdata.c_str());


	//outfile<<"\t" + Utility::Join(rsids, "\t")<<"\n";
	std::vector<Individual>::iterator itr = dataset.begin();
	std::vector<Individual>::iterator end = dataset.end();
	outfile<<"RSID\t" + Utility::Join(dataset, "\t")<<"\n";

	uint rcount = rsids.size();
	for (uint i=0; i<rcount; i++) {
		outfile<<rsids[i]<<"\t";
		itr = dataset.begin();
		while (itr != end) 
			outfile<<itr++->genotypes[i]<<" ";
		outfile<<"\n";
	}
	
	outfile.close();
	std::cerr<<"Genotypes written to file: "<<snpdata<<"\n";

	std::string mapFile = std::string(filename) + ".markers";
	outfile.open(mapFile.c_str());

	std::vector<RsID>::iterator ritr = rsids.begin();
	std::vector<RsID>::iterator rend = rsids.end();
	while (ritr != rend) {
		outfile<<ritr->chrom<<" "<<ritr->pos<<" "<<ritr->rsid<<" "<<ritr->ref<<" "<<ritr->alt<<"\n";
		ritr++;
	}
	outfile.close();
	std::cerr<<"Marker info written to: "<<mapFile<<"\n";
	std::cerr<<"A total of "<<rsids.size()<<" markers were found.\n";
}

int main(int argc, char **argv) {
	ChromosomeContainer chromosomes(30);

	if (argc > 1) {
		ParseRegionDefinition(chromosomes, argv[1]);

		ParseVCF(chromosomes, argv[2]);
	} else {
		std::cerr<<"vcf_selective_parser v1.0 -- Report genotypes from vcf file that are contained" <<
			 " within the segments provided by the segment file.\n" <<
			 "Args: \n" <<
			 "\tsegment_filename - This is the name of the file containing the segments of interest\n" <<
			 "\tvcf_filename - The file containing the genotypes\n";
		return 1;
	}

	return 0;
}
