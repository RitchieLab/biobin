/* 
 * File:   ldspline.h
 * Author: torstees
 *
 * Created on August 30, 2010, 12:16 PM
 */

#ifndef LDSPLINE_H
#define	LDSPLINE_H

#include <map>
#include <vector>
#include <fstream>
#include "utility/strings.h"
#include "locuslookup.h"

namespace Spline {

class LdSpline {
public:
	LdSpline();
	LdSpline(const LdSpline& orig);
	virtual ~LdSpline();

	/**
	 * @brief Scan hapmap files for data point counts so that we can build our binary structures
    * @param chrom
    * @param filename
    */
	void LoadHeadersFromHapmap(const char *chrom, const char* filename);

	/**
	 * @brief Intiates the actual LD loading
    * @param filename
    */
	void LoadLdFromHapmap(const char *filename);

	/**
	 * Intiates complete loading from binary. This is probably not what you want, unless you have tons of RAM
    * @param filename
    */
	void LoadFromBinary(const char *filename);

	/**
	 * Writes splines to binary file
    * @param filename
    */
	void SaveToBinary(const char *filename);

	/**
	 * Saves the modifications to disk under a new filename. This really should be a
	 * new filename-because the old file will probably still hold a lot of data.
    * @param newFilename
    */
	void SaveToCopyBinary(const char *newFilename);

	/**
	 * Generate a basic report to the stream
    * @param os
    */
	void RunReport(std::ostream& os);

	/**
	 * Initiate a selective report based on chromosome, position, ld type and value to stream
    * @param chrom
    * @param position
    * @param value
    * @param ldType
    * @param os
    */
	void RunReport(const char *chrom, int position, float value, const char *ldType, std::ostream& os);

	/**
	 * For the spline associated with position, this reports on the min/max loci associated with it
    * @param chrom
    * @param position
    * @param value
    * @param type
    * @param os
    */
	void ReportBounds(const char *chrom, int position, float value, const char *type, std::ostream& os);

	/**
	 * Initiates the haploview conversion using the chromosome/filename mapping found at filename
    * @param filename Format: chromsome hapmap-filename
    * @return Name of the new binary representation
    */
	std::string ConvertHaploLD(const char *filename);

	/**
	 * Exports BIM file for use with Life Over
    * @param bimOrig file name to be written to
    */
	void ExportForLiftOver(const char *bimOrig);

	/**
	 * Imports BIM files containing new BP locations after being updated by life over
    * @param loFilename The file containing new locations
    * @param loUnmapped The unmapped file
    */
	void ImportLiftOver(const char *loFilename, const char *loUnmapped);

	/**
	 * Opens up the binary file (this will not load anything except the most basic header details. The file
	 * should remain open during the use since things are loaded on demand.
    * @param filename
    * @param loadFullHeaders If true, the entire header will be loaded (which is usually not needed,
	 * and requires extra RAM and time)
    */
	void OpenBinary(const char *filename, bool loadFullHeaders = false);

	/**
	 * Generates a quick summary of the contents at chromosome, chrom
    * @param os
    * @param chrom
    */
	void Summarize(std::ostream& os, const char *chrom);

	/**
	 * Returns 0 or more splines that fall within the bp range start and stop. These
	 * two values do not have to exist as loci within the spline file
    * @param chrom
    * @param start Lower bound position for the splines
    * @param stop Upper bound positions for the splines
    * @return
    */
	std::vector<SnpSpline> GetLocusRange(const char *chrom, int start, int stop);

	/**
	 * Returns the spline boundaries for a given DP value. Upper and lower bounds don't have to
	 * exist in the file for this function to work.
    * @param chrom
    * @param lowerBound Lower bound bp location for the search
    * @param upperbound Upper bound bp location for the search
    * @param value
    * @return
    */
	std::pair<int, int> GetBoundariesDP(const char *chrom, int lowerBound, int upperbound, float value);

	/**
	 * Returns the spline boundaries for a given RS value
    * @param chrom
    * @param lowerBound Lower bound bp location for the search
    * @param upperbound Upper bound bp location for the search
    * @param value
    * @return
    */
	std::pair<int, int> GetBoundariesRS(const char *chrom, int lowerBound, int upperbound, float value);

	/**
	 * Returns the spline lookup for a given chromosome. This is useful for repeated queries
    * @return
    */
	std::map<std::string, LocusLookup> GetChromosomes();
private:
	std::map<std::string, LocusLookup> loci;										///< This is how we record our loci. The spline will return only indexes
	std::map<std::string, std::string> filenames;								///< Used to record the hapmap files so that we can parse them for LD after we've set up the headerspace
	std::fstream file;																	///< This is the file that is or will hold the data
};

}
#endif	/* LDSPLINE_H */

