#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include "snpspline.h"
#include <string>
#include <set>

namespace Spline {

class LocusLookup {
public:
	typedef std::map<int, float> SplineLookup;		///< position -> ld statistic

	typedef std::map<int, int> RsToPos;					///< rs -> pos
	typedef std::map<int, int> PosToRS;					///< pos-> rs

	LocusLookup(std::fstream* file = NULL, const char *chromosome = "");
	~LocusLookup();

	/**
	 * Add a locus to the lookup object if the bp location doesn't already exist
    * @param rs rsid of the SNP
    * @param pos bp location
    * @param offset File offset where the spline details will be found
    */
	void AddLocus(int rs, int pos, uint64_t offset);

	/**
	 * Returns the locus index for a given BP location
    * @param pos
    * @return
    */
	int GetIndex(int pos);

	/**
	 * Return the RSID associated with a given index
    * @param idx
    * @return
    */
	int GetRS(int idx);

	/**
	 * Returns the bp location at a given index
    * @param idx
    * @return
    */
	int GetPos(int idx);

	/**
	 * Returns the SNP count associated with the current lookup structure
    * @return
    */
	int GetSnpCount();

	/**
	 * Increments the spline counts between two loci (given by bp location). This
	 * is necessary in order to predict what the offsets will be in the binary file.
    * @param pos1
    * @param pos2
    */
	void IncrementSplineCounts(int pos1, int pos2);

	/**
	 * Associates a DP and RSquared with the two bp locations
    * @param pos1
    * @param pos2
    * @param dp
    * @param rs
    */
	void AddLdValue(int pos1, int pos2, float dp, float rs);

	/**
	 * Dumps the binary details to file.
    * @param offset This is the file offset associated with the beginning of this
	 * set of splines. This is tricky, because each set of splines can
	 * estimate how large the local set will take-but it has to do so within
	 * the context of the previous N splines, of which it has zero information.
    * @return The file offset where the next spline set will be working
    */
	int64_t DumpBinaryHeader(int64_t offset);

	/**
	 * Dumps header details to a new file as a copy
    * @param file
    * @param offset
    * @return
    */
	int64_t DumpBinaryHeader(std::fstream* file, int64_t offset);

	/**
	 * Initiates the actual file dump of the spline information
    */
	void DumpBinary();

	/**
	 * Initiates the actual file dump, but to a new file
    * @param file
    */
	void DumpBinary(std::fstream* file);

	/**
	 * Returns upper and lower bounds (bp locations) for given position and given RS value
    * @param pos
    * @param value
    * @return
    */
	std::pair<int, int> GetLdSplineBoundsRS(int pos, float value);

	/**
	 * Returns upper and lower bounds (bp locations) for given position and given DP value
    * @param pos
    * @param value
    * @return
    */
	std::pair<int, int> GetLdSplineBoundsDP(int pos, float value);

	/**
	 * Returns boundary extensions for a range, where start and stop are BP bounds and might not
	 * actually be found inside the spline file
    * @param start
    * @param stop
    * @param value
    * @return
    */
	std::pair<int, int> GetRangeBoundariesRS(int start, int stop, float value);
	/**
	 * Returns boundary extensions for a range, where start and stop are BP bounds and might not
	 * actually be found inside the spline file
    * @param start
    * @param stop
    * @param value
    * @return
    */
	std::pair<int, int> GetRangeBoundariesDP(int start, int stop, float value);

	/**
	 * Loads the binary header details
    */
	void LoadBinaryHeader();
	void LoadBinary();

	/**
	 * Scans hapmap files and builds spline details (this is a first pass, and will not
	 * capture any LD values). This first pass IS REQUIRED
    * @param filename
    */
	void LoadHeadersFromHapmap(const std::string& filename);

	/**
	 * Load actual pairwise LD from the hapmap files
    * @param filename
    */
	void LoadLdFromHapmap(const std::string& filename);
	/**
	 * Writes out locus details for lift over in BIM format
    * @param os
	 * chr bp-start bp-stop
    */
	void WriteMapForLiftOver(std::ostream& os);

	/**
	 * Parses file containg liftover BIM data, pausing file read when the original location is
	 * found in dropped positions. In cases where a position was dropped from the new build,
	 * position is set to -1
    * @param infile
    * @param droppedPositions
    */
	void LoadMapFromLiftOver(std::istream& infile, std::set<int>& droppedPositions);

	/**
	 * indicate that the client is no longer using the lookup object.
    */
	void Release();

	/**
	 * Performs a bare-bones load of the binary headers. This is necessary to save
	 * time and minimize memory footprint
    */
	void SkimBinaryHeader();

	/**
	 * Loads spline data (completes the header loading, basically). This is necessary
	 * for the lookup to do it's job
    */
	void LoadLocusDetails();

	/**
	 * Returns the chromosome as it is named
    * @return
    */
	std::string Chromosome() const;

	/**
	 * Summarize the contents of the spline, providing details such as chromosome, rs/bp location
	 * and counts of spline details, among other things
    * @param os
    */
	void Summarize(std::ostream& os);

	/**
	 * Returns Splines for the range between start and stop. Start/stop need not be present themselves
    * @param start
    * @param stop
    * @return
    */
	std::vector<SnpSpline> GetLocusRange(int start, int stop);
	/**
	 * Returns the spline range in positions for dprime, dp
    * @param pos - the position of the item you are
    * @return
    */
	SplineLookup GetLdSplineDP(int pos, float dp);
	/**
	 * Returns the spline range in positions for rsquared, rs
    * @param pos - the position of the item you are
    * @return
    */
	SplineLookup GetLdSplineRS(int pos, float rs);

	/**
	 * Conveniently converts rs to position
    * @return
    */
	RsToPos GetRsToPos();			///< Returns rs->pos

	/**
	 * Conveniently convers position to rs
    * @param Pos
    * @return
    */
	int GetPosToRS(int Pos);		///< Returns the RS number or -1

	/**
	 * Returns the lookup conversion map associated with position to rs
    * @return
    */
	PosToRS GetPosToRS();			///< Returns pos->rs

	int LocusCount();

protected:
	std::string chromosome;									///< The chromosome name
	std::vector<SnpSpline> loci;							///< The actual set of splines
	std::map<int, int> posToIdx;							///< index lookup pos->idx
	std::fstream* file;										///< The file coontaining all of this information (or the destination file)
	int locusCount;											///< Number of loci
	std::fstream::off_type headerLoadPosition;		///< Offset for header details associated with a header skim
};

}

