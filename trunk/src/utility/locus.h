/* 
 * File:   locus.h
 * Author: torstees
 *
 * Created on July 1, 2011, 3:45 PM
 */

#ifndef LOCUS_H
#define	LOCUS_H

#include <string>
#include <map>
#include "types.h"
#include "strings.h"
#include <boost/regex.hpp>

const boost::regex chrom_ex("[CHROMSE-_ ]([0123456789MTXY]+)");
const std::string chrom_format("\\1");
namespace Utility {
	


/**
 * These are some legacy functions used by biofilter before creating the more
 * sophisticated Locus object. Not enough time to properly replace these and test 
 * them....
 */
	
/**
 * Translate string representation to compact index
 * @param chr
 * @return
 */
char ChromToInt(const char* chr);

std::string ChromFromInt(int chrom);
std::string ChromFromIntChr(int chrom);
/**
 * Convert String rs to an integer ID
 * @param rsid
 * @return
 * 
 * This is a legacy function for certain biofilter functions that I don't have
 * time to refactor, where the underlying data really is an RS number (at least
 * as it is currently designed). This includes stuff like the binary file format
 * as well as the hapmap data import for ldspline
 */
uint _RSID(const std::string& rsid);
/**
 * Convert integer RS number to a string
 * @param rs
 * @return
 */
std::string _RSID(uint rs);
/**
 * Convert string list of rsIDs into a list of RS numbers
 * @param rsids
 * @return
 * Currently this assumes that the contents are either numerical representations of RS or are prefixed with the letters RS/rs
 */
Utility::StringArray RsIDs(const std::string& rsids, const std::string& sep);

/** End legacy stuff */



struct Allele {
	Allele() : allele(""), freq(0.0) {}
	Allele(std::string& allele, float freq) : allele(allele), freq(freq) { }
	~Allele() {}

	std::string allele;
	float freq;
};

inline
bool AlleleGT(const Allele& a1, const Allele& a2) {
	if (a1.freq==a2.freq)
		return a1.allele > a2.allele;
	return a1.freq > a2.freq;
}

struct Locus {
	Locus() : chrom(-1), pos(-1), role(0), rsid("")  {}
	Locus(char chrom, uint pos, const std::string& id, int role = 0) : chrom(chrom), pos(pos), role(role) {
		RSID(id.c_str());
	}
	Locus(const Locus& other) : chrom(other.chrom), pos(other.pos), role(other.role), rsid(other.rsid) {
		alleles = other.alleles;
	}

	void Sort() {
		std::sort(alleles.begin(), alleles.end(), AlleleGT);
	}

	float MajorAlleleFreq() const {
		return alleles[0].freq;
	}

	float MinorAlleleFreq() const {
		return 1.0 - alleles[0].freq;
	}

	bool IsMinor(std::string& allele) {
		return (allele == alleles[0].allele);
	}
	
	void AddAllele(std::string& allele, float freq) {
		alleles.push_back(Allele(allele, freq));
	}
	
	char EncodeGenotype(uint a1, uint a2) {
		// We'll use the form that keeps heterozygotes separate
		return (char)(a1 * alleles.size() + a2);
	}
	
	void Print(std::ostream& os, const char* del) const;
	
	
	bool operator<(const Locus& other) const;
	
	std::string Chrom() const;
	
	static char _GuessChromosome(const std::string& chr) ;
	
	/**
	 * Calculate distance between two SNPs
	 * @param other
	 * @return -1 if on separate chromosomes-otherwise, it's the absolute difference
	 */
	uint Distance(Locus& other);	
	
	bool IsExon();
	
	bool IsIntron();
	
	bool IsRegulatory();
	
	/**
	 * This is a way to bring in either an RS Number (or integer version) and 
	 * convert it to an RS ID. This is only used for backward compatibility with
	 * the old approach to biofilter SNPs
	 * 
	 * **WARNING** This function assumes that anything that starts with a number is
	 * an RS number without the R and S prefix. If this isn't what you mean, don't
	 * use it!
    */
	static std::string _RSID(const char *rsid);
	static const uint MaxChromosomes = 26;
	
	std::string RSID(const char *rsid);
	const std::string& RSID() const { return rsid; }

	std::vector<Allele> alleles;
	char chrom;
	uint pos;
	int role;
	
	static Utility::IdCollection ExonRoles;
	static Utility::IdCollection IntronRoles;
	static Utility::IdCollection RegulatoryRoles;
protected:
	std::string rsid;
};

typedef std::vector<Locus> SnpArray;

inline
char Locus::_GuessChromosome(const std::string& chr) {
	std::string chrom				= Utility::UpperCase(chr.c_str());
	std::string c					= boost::regex_replace(chrom, chrom_ex, chrom_format, boost::match_default|boost::format_sed);
	return ChromToInt(c.c_str());
}

inline
bool Locus::IsExon() {
	return ExonRoles.find(role) != ExonRoles.end();
}

inline
bool Locus::IsIntron() {
	return IntronRoles.find(role) != IntronRoles.end();
}

inline
bool Locus::IsRegulatory() {
	return RegulatoryRoles.find(role) != RegulatoryRoles.end();
}



inline
Utility::StringArray RsIDs(const std::string& rsnumbers, const std::string& sep) {
	StringArray ids = Split(rsnumbers.c_str(), sep.c_str());
	StringArray rsids;

	StringArray::iterator itr = ids.begin();
	StringArray::iterator end = ids.end();

	while (itr != end) {
		rsids.push_back(Locus::_RSID(itr++->c_str()));
	}
	return rsids;
}


inline
std::string ChromFromIntChr(int chrom) {
	static std::string labels[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X ", "Y ", "XY", "MT"};
	return std::string("chr") + labels[chrom];
}

inline
std::string ChromFromInt(int chrom) {
	static std::string labels[] = {"1 ", "2 ", "3 ", "4 ", "5 ", "6 ", "7 ", "8 ", "9 ", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X ", "Y ", "XY", "MT"};
	return labels[chrom];
}

inline
char ChromToInt(const char* chrom) {
	std::string chr(chrom);
	if (chrom[0] == 'C' || chrom[0] == 'c') {
		chr = chr.erase(0, 3);
	}
	int c = atoi(chr.c_str());
	if (c==0) {
		std::string chrom = boost::to_upper_copy(std::string(chr.c_str()));
		if (chrom.find("X") != std::string::npos)
			c = 23;
		else if (chrom.find("Y") != std::string::npos)
			c = 24;
		else if (chrom.find("XY") != std::string::npos || chrom.find("X|Y") != std::string::npos)
			c = 25;
		else if (chrom == "MT" || chrom == "M")
			c = 26;
	}
	return (char)(c);
}

/**
 * The following two functions shouldn't be used unless absolutely necessary
 */
inline
std::string _RSID(uint rs) {
	std::stringstream ss;
	ss<<"RS"<<rs;
	return ss.str();
}
inline
uint _RSID(const std::string& rsid) {
	std::string rs = rsid;
	if (rs.find("r") != std::string::npos || rs.find("R") != std::string::npos)
			rs.erase(0,2);
	return atoi(rs.c_str());

}







inline
std::string Locus::_RSID(const char *rs) {
	std::string rsid(rs);
	if (rs[0] >= '0' && rs[0] <= '9')
		rsid = std::string("RS")+rsid;
	else
		rsid =  Utility::UpperCase(rs);
	return rsid;
}

inline
std::string Locus::RSID(const char *rs) {
	rsid =  Utility::UpperCase(rs);
	return rsid;
}

inline
bool Locus::operator<(const Locus& other) const {
	if (chrom == other.chrom)
		return pos < other.pos;
	else
		return chrom < other.chrom;
}

inline
uint Locus::Distance(Locus& other) {
	if (chrom == other.chrom)
		return abs((int)(pos-other.pos));

	return (uint)-1;
}

inline
void Locus::Print(std::ostream& os, const char* del) const {
	os<<Chrom()<<del
			  <<pos<<del
			  <<rsid<<del;
	std::vector<Allele>::const_iterator itr = alleles.begin();
	std::vector<Allele>::const_iterator end = alleles.end();
	while (itr != end) {
		os<<itr->allele<<":"<<itr->freq<<del;
		itr++;
	}
}

inline
std::string Locus::Chrom() const {
	return ChromFromIntChr(chrom - 1);
}

} //Utility

#endif	/* LOCUS_H */

