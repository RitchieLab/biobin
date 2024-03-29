/*
 * vcf_file.h
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */

#ifndef VCF_FILE_H_
#define VCF_FILE_H_

#include <bitset>
#include <cmath>
#include <cstdio>
#include <deque>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "output_log.h"
#include "parameters.h"
#include "vcf_entry.h"

using namespace std;

/**
 * Ritchie labs - I've moved these into a vcf namespace to clean up global - E. Torstenson
 */

namespace VCF {
class vcf_file
{
public:
	vcf_file(const string &filename, bool compressed=false, const string &chr="", const string &exclude_chr="");
	~vcf_file();

	const string filename;
	bool compressed;
	vector<string> meta;
	vector<string> indv;
	unsigned int N_indv;

	unsigned int N_entries;

	deque<streampos> entry_file_locations;
	void get_vcf_entry(unsigned int entry_num, string &out);

	vector<bool> include_indv;
	vector<bool> include_entry;

	//This is a change - we are NOT including the genotype on a variant-by-variant basis!
	vector<bool> include_genotype;
	//vector<vector<bool> > include_genotype;

	void apply_filters(const parameters &params);

	void filter_sites(const set<string> &snps_to_keep, const string &snps_to_keep_file, const string &snps_to_exclude_file, bool keep_then_exclude = false);
	void filter_sites_to_keep(const set<string> &snps_to_keep, const string &snps_to_keep_file);
	void filter_sites_to_exclude(const string &snps_to_exclude_file);
	void filter_sites_by_position(const string &chr, int start_pos, int end_pos);
	void filter_sites_by_positions(const string &positions_file);
	void filter_sites_by_quality(double min_quality);
	void filter_sites_by_mean_depth(double min_mean_depth, double max_mean_depth);
	void filter_sites_by_frequency_and_call_rate(double maf, double max_maf, double non_ref_af, double max_non_ref_af, double min_site_call_rate);
	void filter_sites_by_number_of_alleles(int min_alleles, int max_alleles);
	void filter_sites_by_HWE_pvalue(double min_HWE_pvalue);
	void filter_sites_by_BED_file(const string &bed_file, bool BED_exclude = false);
	void filter_sites_by_mask(const string &mask_file, bool invert_mask = false, int min_kept_mask_value=0);
	void filter_sites_by_filter_status(const set<string> &filter_flags_to_remove, const set<string> &filter_flags_to_keep, bool remove_all = false);
	void filter_sites_by_phase();

	void filter_individuals(const set<string> &indv_to_keep, const set<string> &indv_to_exclude, const string &indv_to_keep_filename, const string &indv_to_exclude_filename, bool keep_then_exclude=true);
	void filter_individuals_by_keep_list(const set<string> &indv_to_keep, const string &indv_to_keep_filename);
	void filter_individuals_by_exclude_list(const set<string> &indv_to_exclude, const string &indv_to_exclude_filename);
	void filter_individuals_by_call_rate(double min_call_rate);
	void filter_individuals_by_mean_depth(double min_mean_depth, double max_mean_depth);
	void filter_individuals_by_phase();

	//void filter_genotypes_by_quality(double min_genotype_quality);
	//void filter_genotypes_by_depth(int min_depth, int max_depth);
	//void filter_genotypes_by_filter_flag(const set<string> &filter_flags_to_remove, bool remove_all = false);

	void output_frequency(const string &output_file_prefix, bool output_counts=false, bool suppress_allele_output=false);
	void output_individuals_by_mean_depth(const string &output_file_prefix);
	void output_site_depth(const string &output_file_prefix, bool output_mean=true);
	void output_genotype_depth(const string &output_file_prefix);
	void output_het(const string &output_file_prefix);
	void output_hwe(const string &output_file_prefix);
	void output_SNP_density(const string &output_file_prefix, int bin_size);
	void output_missingness(const string &output_file_prefix);
	void output_genotype_r2(const string &output_file_prefix, int snp_window_size, int bp_window_size, double min_r2);
	void output_interchromosomal_genotype_r2(const string &output_file_prefix, double min_r2=0.1);
	void output_haplotype_r2(const string &output_file_prefix, int snp_window_size, int bp_window_size, double min_r2);
	void output_singletons(const string &output_file_prefix);
	void output_TsTv(const string &output_file_prefix, int bin_size);
	void output_per_site_nucleotide_diversity(const string &output_file_prefix);
	void output_windowed_nucleotide_diversity(const string &output_file_prefix, int window_size);
	void output_site_quality(const string &output_file_prefix);
	void output_FILTER_summary(const string &output_file_prefix);
	void output_kept_and_removed_sites(const string &output_file_prefix);
	void output_LROH(const string &output_file_prefix);
	void output_indv_relatedness(const string &output_file_prefix);
	void output_PCA(const string &output_file_prefix);

	void output_as_012_matrix(const string &output_file_prefix);
	void output_as_plink(const string &output_file_prefix);
	void output_as_plink_tped(const string &output_file_prefix);
	void output_BEAGLE_genotype_likelihoods(const string &output_file_prefix);
	void output_as_IMPUTE(const string &output_file_prefix);
	void output_as_LDhat_phased(const string &output_file_prefix, const string &chr);
	void output_as_LDhat_unphased(const string &output_file_prefix, const string &chr);
	void output_LDhat_locs_file(const string &output_file_prefix, const string &chr, unsigned int &n_sites_out);
	void output_FORMAT_information(const string &output_file_prefix, const string &FORMAT_id);

	void output_fst(const string &output_file_prefix, vcf_file &vcf_fst);

	void output_sites_in_files(const string &output_file_prefix, vcf_file &diff_vcf_file);
	void output_indv_in_files(const string &output_file_prefix, vcf_file &diff_vcf_file);
	void output_discordance_by_site(const string &output_file_prefix, vcf_file &diff_vcf_file);
	void output_discordance_matrix(const string &output_file_prefix, vcf_file &diff_vcf_file);
	void output_discordance_by_indv(const string &output_file_prefix, vcf_file &diff_vcf_file);
	void output_switch_error(const string &output_file_prefix, vcf_file &diff_vcf_file);

	void output_INFO_for_each_site(const string &output_file_prefix, const vector<string> &INFO_to_extract);

	void print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO=false);

	int N_kept_individuals() const;
	int N_kept_sites() const;
	unsigned int N_genotypes_included(unsigned int entry_num) const;

private:
	ifstream vcf_in;
	gzFile gzvcf_in;
	void open();
	void close();
	bool feof();
	inline void read_line(string &out);
	void parse_header(const string &line);
	void parse_meta(const string &line);
	streampos get_filepos();
	void set_filepos(streampos &filepos);

	bool has_body;
	bool has_file_format;
	bool has_genotypes;
	bool has_header;
	bool has_meta;

	void scan_file(const string &chr="", const string &exclude_chr="");

	void return_indv_union(vcf_file &file2, map<string, pair< int, int> > &combined_individuals);
	void return_site_union(vcf_file &file2, map<pair<string, int>, pair<int, int> > &out);
};
} // Namespace VCF
#endif /* VCF_FILE_H_ */
