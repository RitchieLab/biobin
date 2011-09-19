/*
 * vcf_entry_setters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"
/**
 * Ritchie labs - I've moved these into a vcf namespace to clean up global - E. Torstenson
 */

namespace VCF {

void vcf_entry::set_CHROM(const string &in)
{
	CHROM = in;
}

void vcf_entry::set_POS(const int in)
{
	POS = in;
}

void vcf_entry::set_ID(const string &in)
{
	ID = in;
}

void vcf_entry::set_REF(const string &in)
{
	REF = in;
}

void vcf_entry::set_ALT(const string &in)
{
	istringstream ss(in);
	string tmpstr;
	ALT.resize(0);
	while(!ss.eof())
	{
		getline(ss, tmpstr, ',');
		add_ALT_allele(tmpstr);
	}
	parsed_ALT = true;
}

void vcf_entry::add_ALT_allele(const string &in)
{
	if (in != ".")
	{
		if (find(ALT.begin(), ALT.end(),in) == ALT.end())
		{
			ALT.push_back(in);
		}
	}
	parsed_ALT = true;
}

void vcf_entry::set_QUAL(const double in)
{
	QUAL = in;
}

void vcf_entry::add_FILTER_entry(const string &in)
{
	if ((in != "PASS") && (in != ".") && (in != "0"))
		if (find(FILTER.begin(), FILTER.end(), in) == FILTER.end())
			FILTER.push_back(in);
	parsed_FILTER = true;
}

void vcf_entry::set_FORMAT(const string &in)
{
	FORMAT.resize(0);
	FORMAT_to_idx.clear();

	if (in.size() > 0)
	{
		istringstream ss(in);
		string tmpstr;

		unsigned int pos=0;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ':');
			add_FORMAT_entry(tmpstr, pos);
			pos++;
		}
	}
	parsed_FORMAT = true;
}

void vcf_entry::add_FORMAT_entry(const string &in, unsigned int pos)
{
	FORMAT.push_back(in);
	FORMAT_to_idx[in] = pos;
}

// The following function reads in a genotype from a '0/1'-like string.
// Should handle haploid types to, but NOT polyploidy.
void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const string &in)
{
	size_t pos;
	pos = in.find_first_of("/|");
	ploidy.resize(N_indv);

	string allele1 = in.substr(0,pos);
	string allele2;
	if (pos != string::npos)
	{	// autosome
		ploidy[indv] = 2;
		allele2 = in.substr(pos+1);
		set_indv_PHASE(indv, in[pos]);
	}
	else
	{	// Male chrX, or chrY
		ploidy[indv] = 1;
		allele2 = ".";
		set_indv_PHASE(indv, '|');
	}
	set_indv_GENOTYPE_alleles(indv, make_pair(allele1, allele2));
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase)
{
	set_indv_GENOTYPE_ids(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase)
{
	set_indv_GENOTYPE_alleles(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, const pair<string, string> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (in.first != ".")
		a.first = str2int(in.first);

	if (in.second != ".")
		a.second = str2int(in.second);

	GENOTYPE[indv] = a;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));
	GENOTYPE[indv] = in;
}

void vcf_entry::set_indv_PHASE(unsigned int indv, char in)
{
	if (PHASE.size() == 0)
		PHASE.resize(N_indv, '/');

	PHASE[indv] = in;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GQUALITY(unsigned int indv, double in)
{
	if (GQUALITY.size() == 0)
		GQUALITY.resize(N_indv, -1);

	if (in > 99)
		in = 99;
	GQUALITY[indv] = in;
	parsed_GQ[indv] = true;
}

void vcf_entry::set_indv_DEPTH(unsigned int indv, int in)
{
	if (DEPTH.size() == 0)
		DEPTH.resize(N_indv, -1);

	DEPTH[indv] = in;
	parsed_DP[indv] = true;
}

void vcf_entry::add_indv_GFILTER(unsigned int indv, const string &in)
{
	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	if ((in != ".") && (in != "0"))
		if (find(GFILTER[indv].begin(), GFILTER[indv].end(), in) == GFILTER[indv].end())
			GFILTER[indv].push_back(in);
	parsed_FT[indv] = true;
}

void vcf_entry::set_indv_GFILTER(unsigned int indv, const string &in)
{
	if (in.size() == 0)
		return;
	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	GFILTER[indv].resize(0);

	static istringstream ss;
	static string ith_FILTER;
	ss.clear();
	ss.str(in);
	while (!ss.eof())
	{
		getline(ss, ith_FILTER, ';');

		if ((ith_FILTER == ".") || (ith_FILTER == ""))
			continue;	// Don't bother storing "unfiltered" state.

		GFILTER[indv].push_back(ith_FILTER);
	}
	parsed_FT[indv] = true;
}

void vcf_entry::set_FILTER(const string &FILTER_str)
{
	FILTER.resize(0);
	passed_filters = false;
	if (FILTER_str == "PASS")
		passed_filters = true;
	else
	{
		if ((FILTER_str != "0") && (FILTER_str != "."))
		{
			istringstream ss(FILTER_str);
			string ith_FILTER;
			while (!ss.eof())
			{
				getline(ss, ith_FILTER, ';');
				FILTER.push_back(ith_FILTER);
			}
		}
	}
	parsed_FILTER = true;
}

void vcf_entry::set_INFO(const string &INFO_str)
{
	INFO.resize(0);
	if ((INFO_str.size() > 0) && (INFO_str != "."))
	{
		istringstream ss(INFO_str);
		string tmpstr;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ';');

			istringstream ss2(tmpstr);
			getline(ss2, tmpstr, '=');
			pair<string, string> INFO_entry(tmpstr, ".");

			if (!ss2.eof())
			{	// If there is a value entry, read it now
				getline(ss2, tmpstr);
				INFO_entry.second = tmpstr;
			}
			else	// Otherwise, set it equal to 1
				INFO_entry.second = "1";

			INFO.push_back(INFO_entry);
		}
	}
	parsed_INFO = true;
}

void vcf_entry::add_INFO_descriptor(const string &in)
{
	size_t found=in.find("##INFO=");
	if (found!=string::npos)
	{	// Found an INFO descriptor
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		Field_description I;

		vector<string> tokens;
		tokenize(details, ',', tokens);

		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else error("Expected ID entry in INFO description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Number") I.N_entries =  str2int(entry[1]);
		else error("Expected Number entry in INFO description: " + in);

		tokenize(tokens[2], '=', entry);
		if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") I.Type = Integer;
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) I.Type = Float;
			else if (entry[1] == "Character") I.Type = Character;
			else if (entry[1] == "String") I.Type = String;
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				if (I.N_entries != 0) error("Flag Type must have 0 entries: " + in);
			}
			else error("Unknown Type in INFO meta-information: " + in);
		}
		else error("Expected Type entry in INFO description: " + in);

		tokenize(tokens[3], '=', entry);
		if (entry[0] == "Description")
		{
			I.Description = entry[1];
			for (unsigned int i=4; i<tokens.size(); i++)
			{
				I.Description += "; " + tokens[i];
			}
		}
		else error("Expected Description entry in INFO description: " + in);

		INFO_map[I.ID] = I;
	}
}

void vcf_entry::add_FILTER_descriptor(const string &in)
{
	size_t found=in.find("##FILTER=");
	if (found!=string::npos)
	{
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		vector<string> tokens;
		tokenize(details, ',', tokens);

		string ID, Description;
		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") ID = entry[1];
		else error("Expected ID field in FILTER description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Description")
		{
			Description = entry[1];
			for (unsigned int i=2; i<tokens.size(); i++)
			{
				Description += "; " + tokens[i];
			}
		}
		else
			error("Expected Description field in FILTER description: " + in);

		FILTER_map[ID] = Description;
	}
}

void vcf_entry::add_FORMAT_descriptor(const string &in)
{
	size_t found=in.find("##FORMAT=");
	if (found!=string::npos)
	{	// Found an FORMAT descriptor
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		vector<string> tokens;
		tokenize(details, ',', tokens);
		Field_description I;

		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else error("Expected ID entry in FORMAT description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Number") I.N_entries =  str2int(entry[1]);
		else error("Expected Number entry in FORMAT description: " + in);

		tokenize(tokens[2], '=', entry);
		if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") I.Type = Integer;
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) I.Type = Float;
			else if (entry[1] == "Character") I.Type = Character;
			else if (entry[1] == "String") I.Type = String;
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				if (I.N_entries != 0) error("Flag Type must have 0 entries: " + in);
			}
			else error("Unknown Type in FORMAT meta-information: " + in);
		}
		else error("Expected Type entry in FORMAT description: " + in);

		tokenize(tokens[3], '=', entry);
		if (entry[0] == "Description")
		{
			I.Description = entry[1];
			for (unsigned int i=4; i<tokens.size(); i++)
			{
				I.Description += "; " + tokens[i];
			}
		}
		else error("Expected Description entry in FORMAT description: " + in);

		FORMAT_map[I.ID] = I;
	}
}

} // VCF
