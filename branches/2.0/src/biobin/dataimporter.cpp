/* 
 * File:   dataimporter.cpp
 * Author: torstees
 * 
 * Created on June 16, 2011, 4:21 PM
 */

#include <utility>
#include <algorithm>

#include "dataimporter.h"

#include "knowledge/Locus.h"

using std::pair;
using std::sort;

using Knowledge::Locus;

namespace BioBin {

bool DataImporter::CompressedVCF = false;
bool DataImporter::KeepCommonLoci = true;
bool DataImporter::RareCaseControl = false;

/*bool DataImporter::open(const string& filename) {
	vcf = VCF::vcf_file(filename);
	return true;
}*/

/**
 * TODO ... Um...does snpIndex and line number really mean the same thing...gotta check
 */
void DataImporter::parseSNP(Knowledge::Locus& loc, vector<short>& genotype_map_out) {
	std::string line;
	unordered_map<Knowledge::Locus*, int>::const_iterator pos = _locus_position.find(&loc);
	if(pos != _locus_position.end()){
		int snpIndex = (*pos).second;
		vcf.get_vcf_entry(snpIndex, line);
		//std::cout<<"ParseSNP("<<snpIndex<<") : "<<line<<"\n";
		VCF::vcf_entry entry(vcf.N_indv);
		entry.reset(line);

		// Right now, we'll extract only genotypes. We might
		// want to add in quality or depth or FT later
		entry.parse_genotype_entries(true, false, false, false);
		pair<int,int> genotype;
		genotype_map_out.clear();
		genotype_map_out.reserve(vcf.N_indv);
		// If either of the alleles are -1, then the genotype
		// is char(-1) otherwise, it's the count of major alleles
		for (uint i=0; i<vcf.N_indv; i++) {
			genotype_map_out.push_back(0);
			entry.get_indv_GENOTYPE_ids(i, genotype);
			if (genotype.first == -1 || genotype.second == -1){
				genotype_map_out[i] = -1;
			} else {
				genotype_map_out[i] = loc.encodeGenotype(genotype.first, genotype.second);
			}
		}
	}

}

float DataImporter::getMAF(const vector<int>& allele_count, uint nmcc){
	if(nmcc == 0){
		return -1;
	}else{
		vector<int> ordered_allele_count(allele_count.begin(), allele_count.end());
		sort(ordered_allele_count.begin(), ordered_allele_count.end());
		float ret_val = 0;
		if (ordered_allele_count.size() > 1){
			ret_val = (*(++ordered_allele_count.rbegin())) / ((float) nmcc);
		}
		return ret_val;
//		return ordered_allele_count.size() > 1 ?
//				(*(++ordered_allele_count.rbegin())) / ((float) nmcc) : 0;
	}
}

void DataImporter::remapLocus(Knowledge::Locus* orig_loc, Knowledge::Locus* new_loc){
	unordered_map<Knowledge::Locus*, int>::const_iterator curr_pos = _locus_position.find(orig_loc);
	if(curr_pos != _locus_position.end()){
		_locus_position[new_loc] = (*curr_pos).second;
		_locus_position.erase(curr_pos);
	}

}


}

