/* 
 * File:   genecoverage.h
 * Author: torstees
 *
 * This doesn't fit the same approach as the others, so it shouldn't
 * be planted into the task list. Nonetheless, it is a task.
 * Created on April 14, 2011, 11:58 AM
 */

#ifndef GENE_COVERAGE_H
#define	GENE_COVERAGE_H
#include "task.h"
#include "utility/strings.h"
#include "knowledge/region.h"
#include "liftover/converter.h"

namespace Biofilter {
namespace Task {

	
class GeneCoverage {
public:
	GeneCoverage(Knowledge::SnpDataset* snps, Knowledge::RegionManager *regions, const char *variations);
	GeneCoverage(const GeneCoverage& orig);
	virtual ~GeneCoverage();

	void AddSources(Utility::StringArray& rsFiles, Utility::StringArray& mapFiles, LiftOver::Converter* cnv);
	void GenerateTxtReport(const char *filename);
private:
	void GenerateDetailedTxtReport(const char *filename);
	Knowledge::RegionManager *regions;
	std::string variations;
	Knowledge::SnpDataset *snps;
	std::map<std::string, Knowledge::SnpDataset> datasets;
};



inline GeneCoverage::GeneCoverage(Knowledge::SnpDataset *snps, Knowledge::RegionManager *regions,
		  const char *vars)
	: regions(regions),
		  variations(vars),
		  snps(snps) {}

inline GeneCoverage::GeneCoverage(const GeneCoverage& orig) : regions(orig.regions), variations(orig.variations), snps(orig.snps) { }

inline GeneCoverage::~GeneCoverage() { }


}
}

#endif	/* GENECOVERAGE_H */

