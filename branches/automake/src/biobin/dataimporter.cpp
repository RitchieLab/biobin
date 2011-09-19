/* 
 * File:   dataimporter.cpp
 * Author: torstees
 * 
 * Created on June 16, 2011, 4:21 PM
 */

#include "dataimporter.h"


/***
 * This is going to require some additional work to allow it to be 
 * broken into pieces. Right now, we have no idea which index SNP 0 is
 * for a given chromosome. 
 * 
 * So, for now, this is only going to work if you have a single dataimporter
 * object. However, if you want to use parallel stuff, or deal with different
 * files, that needs to be fixed...
 */

namespace BioBin {

bool DataImporter::CompressedVCF				= false;

/**
 * TODO ... Um...does snpIndex and line number really mean the same thing...gotta check
 */
void DataImporter::ParseSNP(uint snpIndex, std::vector<char>& genotypes) {
	std::string line;
	vcf->get_vcf_entry(snpIndex, line);
	//std::cout<<"ParseSNP("<<snpIndex<<") : "<<line<<"\n";
	entry->reset(line);
	
	// Right now, we'll extract only genotypes. We might 
	// want to add in quality or depth or FT later
	entry->parse_genotype_entries(true, false, false, false);
	pair<int,int> genotype;
	// If either of the alleles are -1, then the genotype
	// is char(-1) otherwise, it's the count of major alleles
	for (uint i=0; i<totalIndividualEntries; i++) {
		genotypes[i] = 0;
		entry->get_indv_GENOTYPE_ids(i, genotype);
		if (genotype.first == -1 || genotype.second == -1)
			genotypes[i] = (char)-1;
		else {
			genotypes[i] = loci[snpIndex].EncodeGenotype(genotype.first, genotype.second);
		}
	}

}

void DataImporter::SetChromosome(char chrom) {
	Close();
	Open(chrom);
}

void DataImporter::GetAllAlleleFrequencies(std::vector<Utility::Locus>& loci) {
	char cachedChrom	= chromosome;
	Open(-1);
	
	std::set<std::string> unknownChromosomes;
	
	//std::vector<Utility::Locus> loci;			///<  Local loci for a given set of 
	uint totalSiteCount					= vcf->N_entries;
	std::string line;
	std::vector<int> alleleCounts;
	double nonMissingChrCount			= 0.0;
	uint nmcc								= 0;					///< Just to avoid redundant conversions
	for (uint i=0; i<totalSiteCount; i++) {
		vcf->get_vcf_entry(i, line);
		entry->reset(line);
		entry->parse_basic_entry(true);
		entry->parse_genotype_entries(true);
		uint alleleCount					= entry->get_N_alleles();
		entry->get_allele_counts(alleleCounts, nmcc, vcf->include_indv, vcf->include_genotype[i]);
		nonMissingChrCount				= nmcc;
		int pos								= entry->get_POS();
		std::string	rsid					= entry->get_ID();
		char chr								= Utility::Locus::_GuessChromosome(entry->get_CHROM());
		Utility::Locus loc				= Utility::Locus(chr, pos, rsid);

		if (chr == 0) // TODO Determine how to handle these that we don't recognize. We need to avoid pulling them when we pull genotypes
			unknownChromosomes.insert(entry->get_CHROM());

		std::string allele				= entry->get_REF();
		double freq							= alleleCounts[0] / nonMissingChrCount;
		loc.AddAllele(allele, freq);
		
		//From here, they are all "ALT" alleles. We would have to evaluate an 
		//if should we want to roll them all in together, since it's a different 
		//call (with a different index scheme...)
		for (uint n = 1; n<alleleCount; n++)	{
			freq								= alleleCounts[n] / nonMissingChrCount;
			allele							= entry->get_ALT_allele(n-1);
			loc.AddAllele(allele, freq);
		}
		
		loc.Sort();
		loci.push_back(loc);
	}
	this->loci.insert(this->loci.end(), loci.begin(), loci.end());	

	Open(cachedChrom);
}

std::vector<Utility::Locus> DataImporter::GetAlleleFrequencies() {
	std::vector<Utility::Locus> loci;			///<  Local loci for a given set of 
	uint totalSiteCount					= vcf->N_entries;
	std::string line;
	std::vector<int> alleleCounts;
	double nonMissingChrCount			= 0.0;
	uint nmcc								= 0;					///< Just to avoid redundant conversions
	for (uint i=0; i<totalSiteCount; i++) {
		vcf->get_vcf_entry(i, line);
		entry->reset(line);
		entry->parse_basic_entry(true);
		entry->parse_genotype_entries(true);
		uint alleleCount					= entry->get_N_alleles();
		entry->get_allele_counts(alleleCounts, nmcc, vcf->include_indv, vcf->include_genotype[i]);
		nonMissingChrCount				= nmcc;
		int pos								= entry->get_POS();
		std::string	rsid					= entry->get_ID();
		Utility::Locus loc				= Utility::Locus(chromosome+1, pos, rsid);

		std::string allele				= entry->get_REF();
		double freq							= alleleCounts[0] / nonMissingChrCount;
		loc.AddAllele(allele, freq);
		
		//From here, they are all "ALT" alleles. We would have to evaluate an 
		//if should we want to roll them all in together, since it's a different 
		//call (with a different index scheme...)
		for (uint n = 1; n<alleleCount; n++)	{
			freq								= alleleCounts[n] / nonMissingChrCount;
			allele							= entry->get_ALT_allele(n-1);
			loc.AddAllele(allele, freq);
		}
		
		loc.Sort();
		loci.push_back(loc);
	}
	this->loci.insert(this->loci.end(), loci.begin(), loci.end());
	return loci;
}

}

