/* 
 * File:   taskgenereport.h
 * Author: torstees
 *
 * Created on April 12, 2011, 5:04 PM
 */
#ifndef TASK_COVERAGE_REPORT_H
#define	TASK_COVERAGE_REPORT_H

#include "task.h"
#include <string>
#include "knowledge/regionmanager.h"
#include "knowledge/snpdataset.h"


namespace Biofilter {
namespace Task {

class CoverageReport : public Task {
public:
	CoverageReport();

	CoverageReport(Utility::StringArray &rsFiles,
		Utility::StringArray &mapFiles,
		const char *variations,
		const char *geneFilename,
		const char *filename);

	CoverageReport(const CoverageReport& orig);

	virtual ~CoverageReport();

	void ExecuteTask();
private:
	Utility::StringArray rsFiles;
	Utility::StringArray mapFiles;
	std::string variations;
	std::string geneFilename;
	std::string filename;
};

inline
CoverageReport::CoverageReport() : Task(1),
		variations(""),
		geneFilename(""),
		filename("") {}

inline
CoverageReport::CoverageReport(Utility::StringArray &rsFiles,
		Utility::StringArray &mapFiles,
		const char *variations,
		const char *geneFilename,
		const char *filename)
			: Task(1),
			  variations(""),
			  geneFilename(geneFilename),
			  filename(filename) {}

inline
CoverageReport::CoverageReport(const CoverageReport& orig)
			: Task(1),
			  variations(""),
			  rsFiles(orig.rsFiles),
			  mapFiles(orig.mapFiles),
			  geneFilename(orig.geneFilename),
			  filename(orig.filename) { }

inline CoverageReport::~CoverageReport() { }

inline void CoverageReport::ExecuteTask() {
	std::map<std::string, Knowledge::SnpDataset> datasets;
	std::set<std::string> rsids;

	Utility::StringArray::iterator itr = rsFiles.begin();
	Utility::StringArray::iterator end = rsFiles.end();

	while (itr != end) {
		std::string filename = *itr;
		std::string fileContents = Utility::LoadContents(filename.c_str());
		Utility::StringArray rsIDs = Utility::RsIDs(fileContents, "\n");
		rsids.insert(rsIDs.begin(), rsIDs.end());
		itr++;
	}

	SnpDataset &master = Knowledge::SnpDataset(variations);

	itr = mapFiles.begin();
	end = mapFiles.end();

	while (itr != end) {
		std::string filename = *itr++;
		abort();
	}

	//OK, now we create the RS datasets
	itr = rsFiles.begin();
	end = rsFiles.end();

	while (itr != end) {
		std::string filename				= *itr;
		std::string fileContents		= Utility::LoadContents(filename);
		Utility::IdCollection rsIDs	= Utility::RsIDs(fileContents, "\n");
		master.GetFragment(rsIDs, datasets[filename]);
		itr++;
	}




}

}

}

#endif	/* TASKGENEREPORT_H */
