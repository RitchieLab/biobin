/* 
 * File:   taskgenereport.h
 * Author: torstees
 *
 * Created on April 12, 2011, 5:04 PM
 */
#ifndef TASKGENEREPORT_H
#define	TASKGENEREPORT_H

#include "taskgenebased.h"
#include <string>


namespace Biofilter {
namespace Task {

class GeneReport : public GeneBased {
public:
	GeneReport();

	virtual ~GeneReport();

	virtual std::string GetFileSuffix() { return "gene-report"; }
	virtual std::string GetFileExtension() { return "csv"; }
	virtual std::string GetFileDescription() { return "Region Details Report";}

	void ExecuteTask();
private:
	std::string filename;
};

inline
GeneReport::GeneReport() : GeneBased() {}


inline
GeneReport::~GeneReport() { }

inline
void GeneReport::ExecuteTask() {
	regions->GenerateGeneReport(filename.c_str(), *snps);
}

}

}

#endif	/* TASKGENEREPORT_H */
