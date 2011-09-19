/* 
 * File:   regionconversion.h
 * Author: torstees
 *
 * Created on April 28, 2011, 2:22 PM
 */

#ifndef REGIONCONVERSION_H
#define	REGIONCONVERSION_H
namespace Knowledge {


class BuildTranslator {
public:




	struct RegionSet {

	};

	BuildTranslator();
	BuildTranslator(const BuildTranslator& orig);
	virtual ~BuildTranslator();
private:
	int id;						///< Chain ID
	long score;					///< Chain Score
	int length;					///< Chain Length
	bool strand;				///< T/F -> +/-
	std::string chrom;		///< Target Chromosome name
};

}
#endif	/* REGIONCONVERSION_H */

