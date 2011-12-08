/* 
 * File:   dataimporter.cpp
 * Author: torstees
 * 
 * Created on June 16, 2011, 4:21 PM
 */

#include <utility>

#include "dataimporter.h"

#include "knowledge/Locus.h"

using std::pair;

/***
 * This is going to require some additional work to allow it to be 
 * broken into pieces. Right now, we have no idea which index SNP 0 is
 * for a given chromosome. 
 * 
 * So, for now, this is only going to work if you have a single dataimporter
 * object. However, if you want to use parallel stuff, or deal with different
 * files, that needs to be fixed...
 */

using Knowledge::Locus;

namespace BioBin {

bool DataImporter::CompressedVCF = false;

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
		genotype_map_out.reserve(vcf.N_indv);
		// If either of the alleles are -1, then the genotype
		// is char(-1) otherwise, it's the count of major alleles
		for (uint i=0; i<vcf.N_indv; i++) {
			genotype_map_out[i] = 0;
			entry.get_indv_GENOTYPE_ids(i, genotype);
			if (genotype.first == -1 || genotype.second == -1){
				genotype_map_out[i] = -1;
			} else {
				genotype_map_out[i] = loc.encodeGenotype(genotype.first, genotype.second);
			}
		}
	}

}



}

