/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */

// Class for reading in, checking and storing user parameters
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <set>

#include "output_log.h"


/**
 * Ritchie labs - I've moved these into a vcf namespace to clean up global - E. Torstenson
 */

namespace VCF {

using namespace std;
	
const string VCFTOOLS_VERSION="v0.1.5";

class parameters
{
public:
	bool BED_exclude;
	string BED_file;
	string chr_to_exclude;
	string chr_to_keep;
	bool diff_discordance_matrix;
	string diff_file;
	bool diff_file_compressed;
	bool diff_indv_discordance;
	bool diff_site_discordance;
	bool diff_switch_error;
	int end_pos;
	string FORMAT_id_to_extract;
	string fst_file;
	bool fst_file_compressed;
	set<string> geno_filter_flags_to_exclude;
	string indv_exclude_file;
	string indv_keep_file;
	set<string> indv_to_exclude;
	set<string> indv_to_keep;
	vector<string> INFO_to_extract;
	set<string> INFO_to_keep;
	bool invert_mask;
	bool keep_all_INFO;
	int ld_bp_window_size;
	int ld_snp_window_size;
	double maf;
	string mask_file;
	int max_alleles;
	int max_genotype_depth;
	double max_indv_mean_depth;
	double max_maf;
	double max_mean_depth;
	double max_non_ref_af;
	int min_alleles;
	int min_genotype_depth;
	double min_genotype_quality;
	double min_HWE_pvalue;
	double min_indv_call_rate;
	double min_indv_mean_depth;
	int min_kept_mask_value;
	double min_mean_depth;
	double min_quality;
	double min_r2;
	double min_site_call_rate;
	double non_ref_af;
	bool output_012_matrix;
	bool output_as_IMPUTE;
	bool output_as_ldhat_phased;
	bool output_as_ldhat_unphased;
	bool output_BEAGLE_genotype_likelihoods;
	bool output_counts;
	bool output_filter_summary;
	bool output_filtered_sites;
	bool output_freq;
	bool output_geno_depth;
	bool output_geno_rsq;
	bool output_hap_rsq;
	bool output_het;
	bool output_HWE;
	bool output_indv_depth;
	bool output_interchromosomal_rsq;
	bool output_LROH;
	bool output_missingness;
	bool output_PCA;
	string output_prefix;
	bool output_relatedness;
	bool output_singletons;
	bool output_site_depth;
	bool output_site_mean_depth;
	bool output_site_pi;
	bool output_site_quality;
	int output_SNP_density_bin_size;
	int output_TsTv_bin_size;
	bool phased_only;
	int pi_window_size;
	bool plink_output;
	bool plink_tped_output;
	string positions_file;
	bool recode;
	bool remove_all_filtered_genotypes;
	bool remove_all_filtered_sites;
	set<string> site_filter_flags_to_exclude;
	set<string> site_filter_flags_to_keep;
	string snps_to_exclude_file;
	string snps_to_keep_file;
	set<string> snps_to_keep;
	int start_pos;
	bool suppress_allele_output;
	string vcf_filename;
	bool vcf_compressed;

	parameters(int argc, char *argv[]);
	~parameters(){};

	void read_parameters();
	void print_help();
	void print_params();

private:
	void check_parameters();
	static void error(string err_msg, int code);

	vector<string> argv;

	string get_arg(unsigned int i);
};

} //VCF
#endif /* PARAMETERS_H_ */
