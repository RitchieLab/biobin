/*
 * File:   taskgenereport.h
 * Author: torstees
 *
 * Created on April 12, 2011, 5:04 PM
 */
#ifndef TASK_SNP_GENE_MAP_H
#define	TASK_SNP_GENE_MAP_H

#include "taskgenebased.h"
#include <string>


namespace Biofilter {
namespace Task {
class SnpGeneMap : public GeneBased {
public:
	SnpGeneMap();

	SnpGeneMap(const SnpGeneMap& orig);

	virtual ~SnpGeneMap();

	virtual std::string GetFileSuffix() { return "snp-gene-map"; }
	virtual std::string GetFileExtension() { return "csv"; }
	virtual std::string GetFileDescription() { return "SNP/Gene Relationship Report";}

	void ExecuteTask();

};

inline
SnpGeneMap::SnpGeneMap() : GeneBased() {}

inline
SnpGeneMap::SnpGeneMap(const SnpGeneMap& orig)
			: GeneBased(orig) { }

inline
SnpGeneMap::~SnpGeneMap() { }

inline
void SnpGeneMap::ExecuteTask() {
	std::ofstream file(filename.c_str());
	file<<"Chrom,RSID,Gene Name,Location w/in Gene\n";
	uint count = snps->Size();
	for (uint i=0; i<count; i++) {
		std::multimap<uint, uint>::iterator itr = geneLookup->lower_bound(i);
		std::multimap<uint, uint>::iterator end = geneLookup->upper_bound(i);

		Utility::Locus &s = (*snps)[i];
		while (itr != end) {
			Knowledge::Region &r = (*regions)[itr++->second];
			file<<s.Chrom()<<","<<s.RSID()<<","<<r.name<<","<<r.DescribeRelationship(s.pos)<<"\n";
		}
	}
}

}
}

#endif	/* TASKGENEREPORT_H */
