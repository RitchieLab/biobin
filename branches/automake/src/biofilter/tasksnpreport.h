/*
 * File:   taskgenereport.h
 * Author: torstees
 *
 * Created on April 12, 2011, 5:04 PM
 */
#ifndef TASK_SNP_REPORT_H
#define	TASK_SNP_REPORT_H

#include "taskgenebased.h"
#include <string>
#include "utility/locus.h"

namespace Biofilter {

namespace Task {

class SnpReport : public GeneBased {
public:
	SnpReport();


	virtual ~SnpReport();

	virtual std::string GetFileSuffix() { return "snp-report"; }
	virtual std::string GetFileExtension() { return "csv"; }
	virtual std::string GetFileDescription() { return "SNP Report";}

	void ExecuteTask();
};

inline
SnpReport::SnpReport() : GeneBased() {}


inline
SnpReport::~SnpReport() { }

inline
void SnpReport::ExecuteTask() {
	std::ofstream file(filename.c_str());
	file<<"Chrom,RSID,Gene Name,(Other Gene Name...)\n";
	uint count = snps->Size();
	for (uint i=0; i<count; i++) {
		std::multimap<uint, uint>::iterator itr = geneLookup->lower_bound(i);
		std::multimap<uint, uint>::iterator end = geneLookup->upper_bound(i);

		Utility::Locus &s = (*snps)[i];
		file<<s.Chrom()<<","<<s.RSID()<<",";

		Utility::StringArray geneNames;
		while (itr != end) {
			Knowledge::Region &r = (*regions)[itr++->second];
			geneNames.push_back(r.name);
		}

		if (geneNames.size() > 0)
			file<<Utility::Join(geneNames, ":");
		file<<"\n";
	}

}
}
}

#endif	/* TASKGENEREPORT_H */
