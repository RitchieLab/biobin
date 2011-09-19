/* 
 * File:   snp.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:28 PM
 */

#ifndef SNP_H
#define	SNP_H

#include <sys/types.h>
#include <vector>
#include "utility/strings.h"
#include "utility/locus.h"

namespace LiftOver {
/*
struct SNP {
	SNP() : chr(0), pos(0), role(0), rsid("") {}
	SNP(char chr, uint pos, const char *rsid, int role = 0) : chr(chr), pos(pos), role(role) {
		//this->rsid = boost::to_upper_copy(std::string(rsid));
		this->rsid = Utility::UpperCase(rsid);
	}
	SNP(const SNP& other) : chr(other.chr), pos(other.pos), role(other.role), rsid(other.rsid) {}

	char chr;				///< Chromosome (integer value 1..25)
	uint pos;				///< position (offset from beginning of chromosome)
	int role;				///< Lookup for the consequence type

	/ **
	 * STL required for sorting
	 * @param other
	 * @return T/F is less than
	 * /
	bool operator<(const SNP& other) const;

	/ **
	 * Calculate distance between two SNPs
	 * @param other
	 * @return -1 if on separate chromosomes-otherwise, it's the absolute difference
	 * /
	uint Distance(SNP& other);

	//std::string RSID(const char *rs) { rsid=boost::to_upper_copy(std::string(rs)); return rsid; }
	std::string RSID(const char *rs) { rsid=Utility::UpperCase(rs); return rsid; }
	std::string RSID() const { return rsid; }
	uint intRSID() const { return Utility::RSID(rsid); }
	
private:
	std::string rsid;		///< Name - this could be almost anything - but we are assuming any letters are upper cased
};
*/
typedef Utility::Locus SNP;
typedef std::vector<Utility::Locus> SnpArray;

/**
inline
bool SNP::operator<(const SNP& other) const {
	if (chr == other.chr)
		return pos < other.pos;
	else
		return chr < other.chr;
}

inline
uint SNP::Distance(SNP& other) {
	if (chr == other.chr)
		return abs((int)(pos-other.pos));

	return (uint)-1;
}
 * */
}
#endif	/* SNP_H */

