#include <map>
#include <vector>
#include <fstream>

namespace Spline {

struct LdStat {
	LdStat(float dp, float rs) : rs(rs), dp(dp) { }
	LdStat() : rs(-1.0), dp(-1.0) { }

	float rs;
	float dp;
};

class SnpSpline {
public:
	SnpSpline(int idx, int rsid, int pos, uint64_t offset);
	~SnpSpline();

	/**
	 * push a new set of DP and RSquared upstream
    */
	void AddStatsUpstream(int idx, float dp, float rs);

	/**
	 * push a new set of DP and RSquared downstream
    */
	void AddStatsDownstream(int idx, float dp, float rs);

	/**
	 *  Return the spline maps associated with DP of minDP
    */
	std::map<int, float> GetSplineDP(float minDP);
	/**
	 *  Return the spline maps associated with RS of minRS
    */
	std::map<int, float> GetSplineRS(float minRS);

	/**
	 * Returns the upper and lower bounds for the spline at value, minDP
    */
	std::pair<int, int> GetSplineBoundsRS(float minDP);
	/**
	 * Returns the upper and lower bounds for the spline at value, minRS
    */
	std::pair<int, int> GetSplineBoundsDP(float minDP);

	bool LoadFromBinary(std::fstream* file);
	void WriteToBinary(std::fstream* file);

	void UpdateUpstreamSplineCount(int idx);
	void UpdateDownstreamSplineCount(int idx);
	void Preload(std::fstream *file);
	void ReleaseStats();

	std::vector<LdStat> GetUpstream();
	std::vector<LdStat> GetDownstream();

	int GetSplineCount();
	void Summarize(std::ostream& os, std::fstream* file, const char *chromosome);
	int pos;
	int rs;
	int idx;
	std::fstream::pos_type offset;
	int usIdx;									///< The maximum index, which will be used to determine how large the upstream vector will be
	int dsIdx;									///< The minimum index, this will be used to determine how large the downstream vector will be

protected:
	int pins;									//If we end up using a buffer, we can use this to note how many pins are holding it up-when this reaches zero, release will empty the spline buffer
	std::vector<LdStat> upstream;
	std::vector<LdStat> downstream;
};

}

