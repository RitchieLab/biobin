/*
 * vcf_file.cpp
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_file.h"

/**
 * Ritchie labs - I've moved these into a vcf namespace to clean up global - E. Torstenson
 */

namespace VCF {
vcf_file::vcf_file(const string &filename, bool compressed, const string &chr, const string &exclude_chr) :
	filename(filename),
	compressed(compressed),
	has_body(false),
	has_file_format(false),
	has_genotypes(false),
	has_header(false),
	has_meta(false)
{
	open();
	scan_file(chr, exclude_chr);
}

vcf_file::~vcf_file()
{
	close();
}

// Parse VCF meta information
void vcf_file::parse_meta(const string &line)
{
	has_meta = true;
	meta.push_back(line);
	size_t found=line.find("##fileformat=");
	if (found!=string::npos)
	{
		has_file_format = true;
		found = line.find_first_of("=");
		string version = line.substr(found+1);
		if (version != "VCFv4.0")
			error("VCF version must be v4.0:\nYou are using version " + version);
	}

	found=line.find("##INFO=");
	if (found!=string::npos)
	{	// Found an INFO descriptor
		vcf_entry::add_INFO_descriptor(line);
	}

	found=line.find("##FILTER=");
	if (found!=string::npos)
	{	// Found a FILTER descriptor
		vcf_entry::add_FILTER_descriptor(line);
	}

	found=line.find("##FORMAT=");
	if (found!=string::npos)
	{	// Found a genotype filter descriptor
		vcf_entry::add_FORMAT_descriptor(line);
	}
}

// Parse VCF header, and extract individuals etc.
void vcf_file::parse_header(const string &line)
{
	// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	(FORMAT	NA00001 NA00002 ... )
	if (has_header == true)
		warning("Multiple Header lines.");

	has_header = true;
	istringstream header(line);
	int count = 0;
	string tmp_str;
	unsigned int N_header_indv = 0;
	has_genotypes = false;
	while (!header.eof())
	{
		header >> tmp_str;
		switch (count)
		{
			case 0: if (tmp_str != "#CHROM") warning("First Header entry should be #CHROM: " + tmp_str); break;
			case 1: if (tmp_str != "POS") warning("Second Header entry should be POS: " + tmp_str); break;
			case 2: if (tmp_str != "ID") warning("Third Header entry should be ID: " + tmp_str); break;
			case 3: if (tmp_str != "REF") warning("Fourth Header entry should be REF: " + tmp_str); break;
			case 4: if (tmp_str != "ALT") warning("Fifth Header entry should be ALT: " + tmp_str); break;
			case 5: if (tmp_str != "QUAL") warning("Sixth Header entry should be QUAL: " + tmp_str); break;
			case 6: if (tmp_str != "FILTER") warning("Seventh Header entry should be FILTER: " + tmp_str); break;
			case 7: if (tmp_str != "INFO") warning("Eighth Header entry should be INFO: " + tmp_str); break;
			case 8:
				if (tmp_str != "FORMAT")
					warning("Ninth Header entry should be FORMAT: " + tmp_str);
				else
					has_genotypes = true;
				break;
			default:
			{
				if (count <= 8)
					error("Incorrectly formatted header.");
				indv.push_back(tmp_str);
				N_header_indv++;
			}
		}
		count++;
	}
	N_indv = N_header_indv;

	if ((has_genotypes == true ) && (N_indv == 0))
		warning("FORMAT field without genotypes?");
}


// Read VCF file
void vcf_file::scan_file(const string &chr, const string &exclude_chr)
{
	printLOG("Scanning " + filename + " ... \n");

	bool filter_by_chr = (chr != "");
	bool exclude_by_chr = (exclude_chr != "");
	string line, tmp;
	N_indv = 0;
	unsigned int N_read = 0;
	istringstream ss;
	string last_CHROM = "";
	N_entries=0;
	string CHROM;
	bool finish = false;
	int last_POS = -1;
	int POS;
	streampos filepos;

	while(!feof())
	{
		filepos = get_filepos();
		read_line(line);

		if (line.length() <= 2)
			continue;

		if (line[0] == '#')
		{
			if (line[1] == '#')
			{	// Meta information
				parse_meta(line);
			}
			else
			{	// Must be header information: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	(FORMAT	NA00001 NA00002 ... )
				parse_header(line);
			}
		}
		else
		{	// Must be a data line
			ss.clear(); ss.str(line);
			ss >> CHROM;

			N_read++;

			if ((filter_by_chr == true) && (last_CHROM == chr) && (CHROM != chr))
			{	// Presuming the file to be sorted (it should be), we have already found the chromosome we wanted, so there's no need to continue.
				printLOG("\tCompleted reading required chromosome. Skipping remainder of file.\n");
				finish = true;
				break;
			}

			if (CHROM != last_CHROM)
			{
				printLOG("Currently scanning CHROM: " + CHROM);
				if ((exclude_by_chr == true) && (CHROM == exclude_chr))
					printLOG(" - excluded.");
				printLOG("\n");
				last_CHROM = CHROM;
				last_POS = -1;
			}

			if ((exclude_by_chr == true) && (CHROM == exclude_chr))
				continue;

			if (filter_by_chr == true)
			{	// For speed, only parse the entry if it's needed
				if (CHROM == chr)
				{
					ss >> POS;
					if (POS < last_POS)
						error("VCF file is not sorted at: " + CHROM + ":" + int2str(POS));
					last_POS = POS;
					entry_file_locations.push_back(filepos);
					N_entries++;
				}
			}
			else
			{
				ss >> POS;
				if (POS < last_POS)
					error("VCF file is not sorted at: " + CHROM + ":" + int2str(POS));
				last_POS = POS;
				entry_file_locations.push_back(filepos);
				N_entries++;
			}
		}
		if (finish == true)
			break;
	}

	if (has_file_format == false)
	{
		warning("*********************************************************************************\n");
		warning("Could not find ##fileformat entry. Assuming VCFv4.0, but THIS MAY CAUSE PROBLEMS.\n");
		warning("*********************************************************************************\n");
	}

	printLOG("\tKeeping " + int2str(N_entries) + " entries (out of " + int2str(N_read) + " read)\n");

	include_indv.clear();
	include_indv.resize(N_indv, true);
	include_entry.clear();
	include_entry.resize(N_entries, true);
	include_genotype.clear();
	include_genotype.resize(N_entries, vector<bool>(N_indv, true));

	printLOG("Done\n");
}

void vcf_file::print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	printLOG("Outputting VCF file... ");
	unsigned int ui;

	string output_file = output_file_prefix + ".recode.vcf";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		error("Could not open VCF Output File: " + output_file, 3);

	for (ui=0; ui<meta.size(); ui++)
		out << meta[ui] << endl;

	out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (N_indv > 0)
		out << "\tFORMAT";
	for (ui=0; ui<N_indv; ui++)
		if (include_indv[ui])
			out << "\t" << indv[ui];
	out << endl;

	string vcf_line;
	for (unsigned int s=0; s<N_entries; s++)
		if (include_entry[s] == true)
		{
			get_vcf_entry(s, vcf_line);
			vcf_entry e(N_indv, vcf_line);
			e.parse_basic_entry(true, true, true);
			e.parse_full_entry(true);
			e.parse_genotype_entries(true,true,true,true);
			e.print(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype[s]);
		}

	out.close();
	printLOG("Done\n");
}

// Return the number of individuals that have not been filtered out
int vcf_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int vcf_file::N_kept_sites() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_entry.size(); ui++)
		if (include_entry[ui] == true)
			N_kept++;
	return N_kept;
}

// Count the number of genotypes that have not been filtered out
unsigned int vcf_file::N_genotypes_included(unsigned int entry_num) const
{
	unsigned int count = 0, ui;
	for (ui=0; ui<N_indv; ui++)
		if ((include_indv[ui] == true) && (include_genotype[entry_num][ui] == true))
		{
			count++;
		}

	return count;
}

void vcf_file::open()
{
	if (!compressed)
	{
		vcf_in.open(filename.c_str(), ios::in);
		if (!vcf_in.is_open())
			error("Could not open VCF file: " + filename, 0);
	}
	else
	{
		gzvcf_in = gzopen(filename.c_str(), "r");
		if (gzvcf_in == NULL)
			error("Could not open GZVCF file: " + filename, 0);
	}
}

void vcf_file::close()
{
	if (!compressed)
		vcf_in.close();
	else
		gzclose(gzvcf_in);
}

bool vcf_file::feof()
{
	bool out;
	if (!compressed)
		out = vcf_in.eof();
	else
	{
		out = gzeof(gzvcf_in);	// Returns 1 when EOF has previously been detected reading the given input stream, otherwise zero.
	}
	return out;
}

streampos vcf_file::get_filepos()
{
	if (!compressed)
		return vcf_in.tellg();
	else
	{
		return gztell(gzvcf_in);	// TODO: Type check
	}
}

void vcf_file::set_filepos(streampos &filepos)
{
	if (!compressed)
	{
		vcf_in.clear();
		vcf_in.seekg(filepos, ios::beg);
	}
	else
	{
		gzseek(gzvcf_in, filepos, SEEK_SET);
	}
}

// TODO: Possibly can make some nice performance improvements here. For example, store the last filepos as a static, and only update the filepos if needed.
void vcf_file::get_vcf_entry(unsigned int entry_num, string &out)
{
	streampos filepos = entry_file_locations[entry_num];
	set_filepos(filepos);
	read_line(out);
}

void vcf_file::read_line(string &out)
{
	if (!compressed)
	{
		getline(vcf_in, out);
		//out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
	else
	{
		static unsigned int MAX_LINE_LEN = 1048576;
		// char * gzgets (gzFile file, char *buf, int len);
		static char *buffer = new char [MAX_LINE_LEN];
		gzgets(gzvcf_in, buffer, MAX_LINE_LEN);
		out = buffer;
		//out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
}
}
