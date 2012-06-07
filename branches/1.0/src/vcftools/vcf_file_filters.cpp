/*
 * vcf_file_filters.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 148 $)
 */

#include "vcf_file.h"

/**
 * Ritchie labs - I've moved these into a vcf namespace to clean up global - E. Torstenson
 */

namespace VCF {
	
void vcf_file::apply_filters(const parameters &params)
{
	// Apply all filters in turn.
	filter_individuals(params.indv_to_keep, params.indv_to_exclude, params.indv_keep_file, params.indv_exclude_file);
	filter_sites(params.snps_to_keep, params.snps_to_keep_file, params.snps_to_exclude_file);
	filter_sites_by_filter_status(params.site_filter_flags_to_exclude, params.site_filter_flags_to_keep, params.remove_all_filtered_sites);
	filter_sites_by_position(params.chr_to_keep, params.start_pos, params.end_pos);
	filter_sites_by_positions(params.positions_file);
	filter_sites_by_BED_file(params.BED_file, params.BED_exclude);
	filter_sites_by_number_of_alleles(params.min_alleles, params.max_alleles);
	filter_sites_by_quality(params.min_quality);
	filter_sites_by_mean_depth(params.min_mean_depth, params.max_mean_depth);
	filter_sites_by_mask(params.mask_file, params.invert_mask, params.min_kept_mask_value);
	filter_individuals_by_mean_depth(params.min_indv_mean_depth, params.max_indv_mean_depth);
	if (params.phased_only == true)
	{
		filter_individuals_by_phase();
		filter_sites_by_phase();
	}
	filter_genotypes_by_quality(params.min_genotype_quality);
	filter_genotypes_by_depth(params.min_genotype_depth, params.max_genotype_depth);
	filter_genotypes_by_filter_flag(params.geno_filter_flags_to_exclude, params.remove_all_filtered_genotypes);
	filter_individuals_by_call_rate(params.min_indv_call_rate);
	filter_sites_by_frequency_and_call_rate(params.maf, params.max_maf, params.non_ref_af, params.max_non_ref_af, params.min_site_call_rate);
	filter_sites_by_HWE_pvalue(params.min_HWE_pvalue);
}

void vcf_file::filter_genotypes_by_quality(double min_genotype_quality)
{
	// Filter genotypes by quality
	if ((min_genotype_quality <= 0) || (has_genotypes == false))
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter genotypes by Quality.");

	printLOG("Filtering out Genotypes with Quality less than " + dbl2str(min_genotype_quality,0) + "\n");
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_genotype_entries(false, true);
		e.filter_genotypes_by_quality(include_genotype[s], min_genotype_quality);
	}
}

void vcf_file::filter_genotypes_by_depth(int min_depth, int max_depth)
{
	// Filter genotypes by depth
	if ((min_depth <= 0) && (max_depth == numeric_limits<int>::max()))
		return;
	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter genotypes by Depth.");

	printLOG("Filtering out Genotypes with Depth less than " + dbl2str(min_depth,0) + " and greater than " + dbl2str(max_depth, 0) + "\n");
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_genotype_entries(false, false, true);
		e.filter_genotypes_by_depth(include_genotype[s], min_depth, max_depth);
	}
}

void vcf_file::filter_genotypes_by_filter_flag(const set<string> &filter_flags_to_remove, bool remove_all)
{
	// Filter genotypes by Filter Flags
	if ((remove_all == false) && (filter_flags_to_remove.size() == 0))
		return;
	if (remove_all == true)
		printLOG("Filtering out all genotypes with FILTER flag.\n");
	else
		printLOG("Filtering out genotypes by Filter Status.\n");

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter genotypes by Filter Flag.");

	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_genotype_entries(false, false, false, true);
		e.filter_genotypes_by_filter_status(include_genotype[s], filter_flags_to_remove, remove_all);
	}
}


void vcf_file::filter_individuals(const set<string> &indv_to_keep, const set<string> &indv_to_exclude, const string &indv_to_keep_filename, const string &indv_to_exclude_filename, bool keep_then_exclude)
{
	// Filter individuals by user provided lists
	if (keep_then_exclude)
	{
		filter_individuals_by_keep_list(indv_to_keep, indv_to_keep_filename);
		filter_individuals_by_exclude_list(indv_to_exclude, indv_to_exclude_filename);
	}
	else
	{
		filter_individuals_by_exclude_list(indv_to_exclude, indv_to_exclude_filename);
		filter_individuals_by_keep_list(indv_to_keep, indv_to_keep_filename);
	}
}

void vcf_file::filter_individuals_by_keep_list(const set<string> &indv_to_keep, const string &indv_to_keep_filename)
{
	// Filter individuals by user provided list
	if ((indv_to_keep_filename == "") && (indv_to_keep.size() == 0))
		return;
	printLOG("Keeping individuals in 'keep' list\n");
	set<string> indv_to_keep_copy = indv_to_keep;
	if (indv_to_keep_filename != "")
	{
		ifstream infile(indv_to_keep_filename.c_str());
		if (!infile.is_open())
			error("Could not open Individual file:" + indv_to_keep_filename, 1);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!infile.eof())
		{
			getline(infile, line);
			ss.str(line);
			ss >> tmp_indv;
			indv_to_keep_copy.insert(tmp_indv);
			ss.clear();
		}
		infile.close();
	}

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (indv_to_keep_copy.find(indv[ui]) == indv_to_keep_copy.end())
			include_indv[ui] = false;
	}
}

void vcf_file::filter_individuals_by_exclude_list(const set<string> &indv_to_exclude, const string &indv_to_exclude_filename)
{
	// Filter individuals by user provided list
	if ((indv_to_exclude_filename == "") && (indv_to_exclude.size() == 0))
		return;
	printLOG("Excluding individuals in 'exclude' list\n");
	set<string> indv_to_exclude_copy = indv_to_exclude;
	if (indv_to_exclude_filename != "")
	{
		ifstream infile(indv_to_exclude_filename.c_str());
		if (!infile.is_open())
		{
			error("Could not open Individual file:" + indv_to_exclude_filename, 1);
		}
		string line;
		string tmp_indv;
		stringstream ss;
		while (!infile.eof())
		{
			getline(infile, line);
			ss.str(line);
			ss >> tmp_indv;
			indv_to_exclude_copy.insert(tmp_indv);
			ss.clear();
		}
		infile.close();
	}
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (indv_to_exclude_copy.find(indv[ui]) != indv_to_exclude_copy.end())
			include_indv[ui] = false;
	}
}

void vcf_file::filter_individuals_by_call_rate(double min_call_rate)
{
	// Filter individuals by call rate
	if (min_call_rate <= 0.0)
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter individuals by call rate.");

	printLOG("Filtering individuals by call rate\n");

	unsigned int ui;
	pair<int, int> genotype;
	vector<int> N_sites_included(N_indv, 0);
	vector<int> N_missing(N_indv, 0);
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, true);
				e.get_indv_GENOTYPE_ids(ui, genotype);
				if (genotype.first != -1)
				{
					N_missing[ui]++;
				}
				N_sites_included[ui]++;
			}
		}
	}

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		double call_rate = N_missing[ui] / (double)N_sites_included[ui];
		if (call_rate < min_call_rate)
			include_indv[ui] = false;
	}
}

void vcf_file::filter_individuals_by_mean_depth(double min_mean_depth, double max_mean_depth)
{
	// Filter individuals by mean depth across sites
	if ((min_mean_depth <= 0) && (max_mean_depth == numeric_limits<double>::max()))
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter individuals by mean depth");

	printLOG("Filtering individuals by mean depth\n");
	unsigned int ui;

	vector<int> N_sites_included(N_indv, 0);
	vector<double> depth_sum(N_indv,0.0);
	int depth;
	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, false, false, true);
				depth = e.get_indv_DEPTH(ui);
				if (depth >= 0)
				{
					depth_sum[ui] += depth;
					N_sites_included[ui]++;
				}
			}
		}
	}

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		double mean_depth = depth_sum[ui] / N_sites_included[ui];
		if ((mean_depth < min_mean_depth) || (mean_depth > max_mean_depth))
			include_indv[ui] = false;
	}
}

void vcf_file::filter_individuals_by_phase()
{
	// Filter individuals that are completely unphased.
	// TODO: Alter this to allow for a max/min level of unphased-ness.
	printLOG("Filtering Unphased Individuals\n");

	if (has_genotypes == false)
		error("Require Genotypes in VCF file to filter by Phase.");

	unsigned int ui, s;
	vector<unsigned int> indv_count(N_indv, 0);
	vector<unsigned int> indv_count_unphased(N_indv, 0);
	string vcf_line;
	vcf_entry e(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);

		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, true);

			indv_count[ui]++;
			if (e.get_indv_PHASE(ui) != '|')
				indv_count_unphased[ui]++;
		}
	}

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		if (indv_count_unphased[ui] == indv_count[ui])
		{
			include_indv[ui] = false;
		}
	}
}


void vcf_file::filter_sites(const set<string> &snps_to_keep, const string &snps_to_keep_file, const string &snps_to_exclude_file, bool keep_then_exclude)
{
	// Filter sites by user provided lists
	if (keep_then_exclude)
	{
		filter_sites_to_keep(snps_to_keep, snps_to_keep_file);
		filter_sites_to_exclude(snps_to_exclude_file);
	}
	else
	{
		filter_sites_to_exclude(snps_to_exclude_file);
		filter_sites_to_keep(snps_to_keep, snps_to_keep_file);
	}
}

void vcf_file::filter_sites_to_keep(const set<string> &snps_to_keep, const string &snps_to_keep_file)
{
	// Filter sites by user provided list
	if ((snps_to_keep.size() == 0) && (snps_to_keep_file == ""))
		return;

	set<string> local_snps_to_keep = snps_to_keep;

	printLOG("Keeping sites by user-supplied list\n");

	if (snps_to_keep_file != "")
	{
		ifstream in(snps_to_keep_file.c_str());
		string tmp;
		if (!in.is_open())
		{
			error("Could not open SNPs to Keep file" + snps_to_keep_file, 0);
		}
		while (!in.eof())
		{
			in >> tmp;
			local_snps_to_keep.insert(tmp);
			in.ignore(numeric_limits<streamsize>::max(), '\n');
		}

		in.close();
	}

	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry();
		if (local_snps_to_keep.find(e.get_ID()) == local_snps_to_keep.end())
			include_entry[s] = false;
	}
}

void vcf_file::filter_sites_to_exclude(const string &snps_to_exclude_file)
{
	// Filter sites by user provided list
	if (snps_to_exclude_file == "")
		return;

	printLOG("Excluding sites by user-supplied list\n");

	set<string> snps_to_exclude;
	if (snps_to_exclude_file != "")
	{
		ifstream in(snps_to_exclude_file.c_str());
		string tmp;
		if (!in.is_open())
		{
			error("Could not open SNPs to Exclude file" + snps_to_exclude_file, 0);
		}
		while (!in.eof())
		{
			in >> tmp;
			snps_to_exclude.insert(tmp);
			in.ignore(numeric_limits<streamsize>::max(), '\n');
		}
		in.close();
	}

	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry();
		if (snps_to_exclude.find(e.get_ID()) != snps_to_exclude.end())
			include_entry[s] = false;
	}
}

void vcf_file::filter_sites_by_quality(double min_quality)
{
	// Filter sites by quality
	if (min_quality < 0)
		return;

	printLOG("Filtering sites with Quality less than " + dbl2str(min_quality,0) + "\n");

	unsigned int s;
	string vcf_line;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;
		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry(true);
		string alt_allele = e.get_ALT_allele(0);
		// The QUAL field has different definitions depending on the state of the
		// alternative allele. Here I treat them separately, although in this case
		// it is unnecessary.
		if ((alt_allele == ".") || (alt_allele == ""))
		{	// The case that the alternative allele is unknown
			// QUAL is -10log_10 p(variant)
			if (e.get_QUAL() < min_quality)
				include_entry[s] = false;
		}
		else
		{	// The normal case
			// QUAL is -10log_10 p(no variant)
			if (e.get_QUAL() < min_quality)
				include_entry[s] = false;
		}
	}
}

void vcf_file::filter_sites_by_mean_depth(double min_mean_depth, double max_mean_depth)
{
	// Filter sites by mean depth
	if ((min_mean_depth <= 0) && (max_mean_depth == numeric_limits<double>::max()))
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file in order to filter sites by mean depth");

	printLOG("Filtering sites by mean depth\n");
	int depth;

	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);

		unsigned int N_indv_included = 0;
		double depth_sum = 0.0;
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e.parse_genotype_entry(ui, false, false, true);
				depth = e.get_indv_DEPTH(ui);
				if (depth >= 0)
				{
					depth_sum += depth;
				}
				N_indv_included++;
			}
		}
		double mean_depth = depth_sum / N_indv_included;

		if ((mean_depth < min_mean_depth) || (mean_depth > max_mean_depth))
			include_entry[s] = false;
	}
}

void vcf_file::filter_sites_by_position(const string &chr, int start_pos, int end_pos)
{
	// Filter sites by user provided position range
	if ((chr == "") || ((start_pos == -1) && (end_pos==numeric_limits<int>::max())))
		return;
	printLOG("Filtering sites by chromosome and/or position\n");
	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;
		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry();
		if (e.get_CHROM() == chr)
		{
			if ((e.get_POS() < start_pos) || (e.get_POS() > end_pos))
				include_entry[s] = false;
		}
		else
			include_entry[s] = false;
	}
}

void vcf_file::filter_sites_by_positions(const string &positions_file)
{
	// Filter sites by a user defined file containing a list of positions
	if (positions_file == "")
		return;
	printLOG("Filtering sites by Positions file\n");
	ifstream BED(positions_file.c_str());
	if (!BED.is_open())
		error("Could not open Positions file: " + positions_file);

	string chr;
	int pos1;
	int idx;
	unsigned int N_chr=0;
	map<string,int> chr_to_idx;
	vector< set<int > > lims;
	stringstream ss;
	string line;
	// Skip header
	BED.ignore(numeric_limits<streamsize>::max(), '\n');
	while (!BED.eof())
	{
		getline(BED, line);
		if (line[0] == '#')
			continue;

		ss.clear();
		ss.str(line);
		ss >> chr >> pos1;

		if (chr_to_idx.find(chr) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[chr] = (N_chr-1);
			lims.resize(N_chr);
		}

		idx = chr_to_idx[chr];
		lims[idx].insert(pos1);
	}
	BED.close();

	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;
		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry();
		e.get_CHROM(chr);
		if (chr_to_idx.find(chr) == chr_to_idx.end())
			include_entry[s] = false;
		else
		{
			pos1 = e.get_POS();
			idx = chr_to_idx[chr];
			bool found=false;

			if (lims[idx].find(pos1) != lims[idx].end())
				found = true;

			if (found == false)
				include_entry[s] = false;
		}
	}
}

void vcf_file::filter_sites_by_BED_file(const string &bed_file, bool BED_exclude)
{
	// Filter sites depending on positions in a BED file.
	if (bed_file == "")
		return;
	printLOG("Filtering sites by BED file\n");
	ifstream BED(bed_file.c_str());
	if (!BED.is_open())
		error("Could not open BED file: " + bed_file);

	string chr;
	int pos1, pos2;
	int idx;
	unsigned int N_chr=0;
	map<string,int> chr_to_idx;
	vector< deque<pair<int,int> > > lims;
	BED.ignore(numeric_limits<streamsize>::max(), '\n');	// Ignore header
	while (!BED.eof())
	{
		BED >> chr >> pos1 >> pos2;
		BED.ignore(numeric_limits<streamsize>::max(), '\n');

		if (chr_to_idx.find(chr) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[chr] = (N_chr-1);
			lims.resize(N_chr);
		}

		idx = chr_to_idx[chr];
		lims[idx].push_back(make_pair(pos1,pos2));
	}
	BED.close();

	for (unsigned int ui=0; ui<lims.size(); ui++)
		sort(lims[ui].begin(), lims[ui].end());

	pair<int,int> range;
	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;
		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry();
		e.get_CHROM(chr);
		if (BED_exclude == false)
		{	// Exclude sites not in BED file
			if (chr_to_idx.find(chr) == chr_to_idx.end())
				include_entry[s] = false;
			else
			{
				pos1 = e.get_POS();
				idx = chr_to_idx[chr];
				bool found=false;
				unsigned int max_ui = lims[idx].size();
				for (unsigned int ui=0; ui<max_ui; ui++)
				{	// TODO: No need to start this loop at zero every time...
					if ((pos1 > lims[idx][ui].first) && (pos1 <= lims[idx][ui].second))
					{
						found=true;
						break;
					}
				}
				if (found == false)
					include_entry[s] = false;
			}
		}
		else
		{	// Exclude sites in BED file
			if (chr_to_idx.find(chr) != chr_to_idx.end())
			{
				pos1 = e.get_POS();
				idx = chr_to_idx[chr];
				bool found=false;
				unsigned int max_ui = lims[idx].size();
				for (unsigned int ui=0; ui<max_ui; ui++)
				{	// TODO: No need to start this loop at zero every time...
					if ((pos1 > lims[idx][ui].first) && (pos1 <= lims[idx][ui].second))
					{
						found=true;
						break;
					}
				}
				if (found == true)
					include_entry[s] = false;
			}
		}
	}
}

void vcf_file::filter_sites_by_mask(const string &mask_file, bool invert_mask, int min_kept_mask_value)
{
	// Filter sites on the basis of a fasta-like mask file.
	if (mask_file == "")
		return;
	if (invert_mask == false)
		printLOG("Filtering sites by mask file\n");
	else
		printLOG("Filtering sites by inverted mask file\n");
	ifstream mask(mask_file.c_str());
	if (!mask.is_open())
		error("Could not open mask file: " + mask_file);

	string line;
	string next_chr="", vcf_line;
	unsigned int next_pos = 0;
	unsigned int next_s = 0;

	unsigned int current_pos = 1;
	string current_header = "";
	bool keep;
	while (!mask.eof())
	{
		getline(mask, line);
		line.erase( line.find_last_not_of(" \t") + 1);

		if (line[0] == '>')
		{	// Header
			current_header = line.substr(1, line.find_first_of(" \t")-1);
			current_pos = 1;
			for (unsigned int s=0; s<N_entries; s++)
			{
				if (include_entry[s] == true)
				{
					get_vcf_entry(s, vcf_line);
					vcf_entry e(N_indv, vcf_line);
					e.parse_basic_entry();
					e.get_CHROM(next_chr);
					if (next_chr == current_header)
					{
						next_pos = (unsigned)e.get_POS();
						next_s = s;
						break;
					}
					else
					{
						include_entry[s] = false;
					}
				}
			}
		}
		else
		{
			if ((current_pos + line.size() >= next_pos) && (next_chr == current_header))
			{
				for (unsigned int ui=0; ui<line.size(); ui++)
				{
					if (current_pos + ui == next_pos)
					{
						char mask_base = line[ui]-48;
						keep = (mask_base <= min_kept_mask_value);
						if (invert_mask == true)
							keep = !keep;

						if (keep == false)
						{
							include_entry[next_s] = false;
						}

						next_s += 1;
						for (unsigned int s=next_s; s<N_entries; s++)
						{
							if (include_entry[s] == true)
							{
								get_vcf_entry(s, vcf_line);
								vcf_entry e(N_indv, vcf_line);
								e.parse_basic_entry();
								e.get_CHROM(next_chr);
								next_pos = (unsigned)e.get_POS();
								next_s = s;
								break;
							}
						}
					}
				}
			}
			current_pos += line.size();
		}
	}
	mask.close();

	// Remaining sites aren't covered by mask, so exclude
	for (unsigned int s=next_s; s<N_entries; s++)
	{
		include_entry[s] = false;
	}
}


void vcf_file::filter_sites_by_number_of_alleles(int min_alleles, int max_alleles)
{
	// Filter sites by the number of alleles (e.g. 2 for bi-allelic)
	if ((min_alleles <= 0) && (max_alleles == numeric_limits<int>::max()))
		return;
	printLOG("Filtering sites by number of alleles\n");

	int N_alleles;
	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);
		e.parse_basic_entry(true);
		N_alleles = e.get_N_alleles();
		if ((N_alleles < min_alleles) || (N_alleles > max_alleles))
		{
			include_entry[s] = false;
		}
	}
}

void vcf_file::filter_sites_by_frequency_and_call_rate(double maf, double max_maf, double non_ref_af, double max_non_ref_af, double min_site_call_rate)
{
	// Filter sites so that all alleles are between limits
	if ((maf <= 0.0) && (max_maf >= 1.0) && (min_site_call_rate <= 0) && (non_ref_af <= 0.0) && (max_non_ref_af >= 1.0))
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file to filter by Frequency and/or Call Rate");

	printLOG("Filtering sites by allele frequency and call rate\n");

	unsigned int N_alleles;
	unsigned int N_non_missing_chr;

	string vcf_line;
	vcf_entry e(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);
		e.parse_basic_entry(true);
		e.parse_genotype_entries(true);
		N_alleles = e.get_N_alleles();

		vector<int> allele_counts;
		e.get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);

		vector<double> freq(N_alleles, 0.0);
		unsigned int Alleles_between_limits = 0;
		for (unsigned int ui=0; ui<N_alleles; ui++)
		{
			freq[ui] = allele_counts[ui] / (double)N_non_missing_chr;
			freq[ui] = min(freq[ui], 1.0 - freq[ui]);

			if ((freq[ui] >= maf) && (freq[ui] <= max_maf))
				Alleles_between_limits++;
		}

		if (Alleles_between_limits != N_alleles)
			include_entry[s] = false;

		//unsigned int N_geno_included = e.get_N_chr();
		double call_rate = N_non_missing_chr / double(e.get_N_chr(include_indv, include_genotype[s]));

		if (call_rate < min_site_call_rate)
			include_entry[s] = false;

		double non_ref_freq = 1.0 - freq[0];
		if ((non_ref_freq < non_ref_af) || (non_ref_freq > max_non_ref_af))
			include_entry[s] = false;
	}
}

void vcf_file::filter_sites_by_HWE_pvalue(double min_HWE_pvalue)
{
	// Filter sites by HWE p-value
	if (min_HWE_pvalue <= 0)
		return;

	if (has_genotypes == false)
		error("Require Genotypes in VCF file to filter sites by HWE.");

	// Note this assumes Biallelic SNPs.
	printLOG("Filtering sites by HWE p-value (only including bi-allelic sites)\n");

	unsigned int b11, b12, b22;
	double p;
	unsigned int N_non_missing_chr;
	vector<int> allele_counts;
	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);

		e.parse_basic_entry(true);
		e.parse_genotype_entries(true);

		e.get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		e.get_genotype_counts(include_indv, include_genotype[s], b11, b12, b22);
		/*
		double freq = allele_counts[0] / (double)N_non_missing_chr;
		double tot = b11 + b12 + b22;
		double exp_11 = freq * freq * tot;
		double exp_12 = 2.0 * freq * (1.0-freq) * tot;
		double exp_22 = (1.0-freq) * (1.0-freq) * tot;

		double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
				+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
				+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22;
		*/

		p = vcf_entry::SNPHWE(b12, b11, b22);

		if (p < min_HWE_pvalue)
		{
			include_entry[s] = false;
		}
	}
}

void vcf_file::filter_sites_by_filter_status(const set<string> &filter_flags_to_remove, const set<string> &filter_flags_to_keep, bool remove_all)
{
	// Filter sites by entries in the FILTER field.
	if ((remove_all == false) && (filter_flags_to_remove.size() == 0) && (filter_flags_to_keep.size() == 0))
		return;

	printLOG("Filtering sites by FILTER Status.\n");

	vector<string> FILTERs;
	string vcf_line;
	unsigned int N_to_remove = filter_flags_to_remove.size();
	unsigned int N_to_keep = filter_flags_to_keep.size();
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_vcf_entry(s, vcf_line);
		vcf_entry e(N_indv, vcf_line);

		e.parse_basic_entry(false, true);

		e.get_FILTER_vector(FILTERs);

		if (N_to_keep > 0)
		{
			bool keep = false;
			for (unsigned int ui=0; ui<FILTERs.size(); ui++)
				if (filter_flags_to_keep.find(FILTERs[ui]) != filter_flags_to_keep.end())
				{
					keep = true; break;
				}

			include_entry[s] = keep;
		}

		if (include_entry[s]==false)
			continue;

		if ((remove_all == true) && (FILTERs.size() > 0))
			include_entry[s] = false;
		else if (N_to_remove > 0)
		{
			for (unsigned int ui=0; ui<FILTERs.size(); ui++)
				if (filter_flags_to_remove.find(FILTERs[ui]) != filter_flags_to_remove.end())
					include_entry[s] = false;
		}
	}

}

void vcf_file::filter_sites_by_phase()
{
	// Filter out sites with unphased entries
	// TODO: Alter this to allow for a max/min level of unphased-ness.
	printLOG("Filtering Sites with Unphased Genotypes\n");
	string vcf_line;
	vcf_entry e(N_indv);

	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		unsigned int count = 0;
		unsigned int count_unphased = 0;
		get_vcf_entry(s, vcf_line);
		e.reset(vcf_line);

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e.parse_genotype_entry(ui, true);

			count++;
			if (e.get_indv_PHASE(ui) != '|')
				count_unphased++;
		}

		if (count_unphased > 0)
			include_entry[s] = false;
	}
}

} //VCF

