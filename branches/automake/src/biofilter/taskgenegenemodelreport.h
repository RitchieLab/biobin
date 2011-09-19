/* 
 * File:   taskgenegenemodelreport.h
 * Author: torstees
 *
 * Created on April 15, 2011, 2:35 PM
 */

#ifndef TASKGENEGENEMODELREPORT_H
#define	TASKGENEGENEMODELREPORT_H
#include "task.h"
#include "knowledge/genegenemodelarchive.h"

namespace Biofilter {
	namespace Task {

class GeneGeneModelReport : public Task {
public:
	GeneGeneModelReport();
	virtual ~GeneGeneModelReport() {}
	virtual void Init(Application* app);
	void ExecuteTask();
protected:
	Knowledge::GeneGeneModelArchive *geneGeneModels;
	std::string geneArchiveName;

};

inline
GeneGeneModelReport::GeneGeneModelReport() : Task(4), geneGeneModels(NULL) {}

inline
void GeneGeneModelReport::Init(Application* app) {
	regions			= app->GetRegions();
	//snps				= app->GetDataset();
	geneGeneModels	= app->GetGeneGeneModels();
	//filename			= app->AddReport("gene-gene-model-report", "csv", "Gene/Gene Model Report");
	filename		= app->AddReport("model-archive", "gene-gene", "Gene/Gene Models");
	geneArchiveName		= app->AddReport("model-archive", "genes", "Gene Definition");
}

inline
void GeneGeneModelReport::ExecuteTask() {
	std::map<float, uint> scores;

	std::cerr<<"Writing Gene/Gene Archive\n";
	regions->WriteArchive(geneArchiveName.c_str(), "\t");
	geneGeneModels->WriteToArchive(filename.c_str(), *regions, scores, Knowledge::BinaryArchive);
	//geneGeneModels->SummarizeModelCounts(scores, *regions);
	std::map<float, uint>::iterator itr = scores.begin();
	std::map<float, uint>::iterator end = scores.end();

	std::cerr<<std::setw(10)<<std::right<<"Impl."<<" "<<std::setw(6)<<std::right<<"Model"<<"\n";
	std::cerr<<std::setw(10)<<std::right<<"Index"<<" "<<std::setw(6)<<std::right<<"Count"<<"\n";
	while (itr != end) {
		std::cerr<<std::setw(10)<<std::right<<itr->first<<" "<<std::setw(6)<<std::right<<itr->second<<"\n";
		itr++;
	}

}

}

}

#endif	/* TASKGENEGENEMODELREPORT_H */

