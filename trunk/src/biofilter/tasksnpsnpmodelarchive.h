/* 
 * File:   taskgenegenemodelreport.h
 * Author: torstees
 *
 * Created on April 15, 2011, 2:35 PM
 */

#ifndef TASKSNPSNPMODELREPORT_H
#define	TASKSNPSNPMODELREPORT_H
#include "task.h"
#include "knowledge/genegenemodelarchive.h"

namespace Biofilter {
	namespace Task {

class SnpSnpModelArchive : public Task {
public:
	SnpSnpModelArchive();
	virtual ~SnpSnpModelArchive() {}
	virtual void Init(Application* app);
	void ExecuteTask();
protected:
	Knowledge::GeneGeneModelArchive *geneGeneModels;
	std::string geneArchiveName;
	Knowledge::SnpSnpModel::Collection snpBasedModels;

};

inline
SnpSnpModelArchive::SnpSnpModelArchive() : Task(4), geneGeneModels(NULL) {}

inline
void SnpSnpModelArchive::Init(Application* app) {
	regions				= app->GetRegions();
	snps					= app->GetDataset();
	geneGeneModels		= app->GetGeneGeneModels();
	//filename			= app->AddReport("gene-gene-model-report", "csv", "Gene/Gene Model Report");
	filename				= app->AddReport("model-archive", "snp-snp", "SNP/SNP Models");
	//geneArchiveName	= app->AddReport("model-archive", "genes", "Gene Definition");
}

inline
void SnpSnpModelArchive::ExecuteTask() {
	std::map<float, uint> scores;

	std::cerr<<"Writing SNP/SNP Model Archive\n";
	geneGeneModels->GenerateModels(snpBasedModels, *regions);
	Knowledge::SnpSnpModel::Collection::iterator itr = snpBasedModels.begin();
	Knowledge::SnpSnpModel::Collection::iterator end = snpBasedModels.end();

	std::ofstream file((filename).c_str());
	if (Knowledge::BinaryArchive) {
		file.close();
		file.open(filename.c_str(), std::ios::binary);
	}
	while (itr != end) {
		scores[itr->ImplicationIndex()]++;
		itr++->Write(file);
	}

	std::map<float, uint>::iterator sitr = scores.begin();
	std::map<float, uint>::iterator send = scores.end();

	std::cerr<<std::setw(10)<<std::right<<"Impl."<<" "<<std::setw(6)<<std::right<<"Model"<<"\n";
	std::cerr<<std::setw(10)<<std::right<<"Index"<<" "<<std::setw(6)<<std::right<<"Count"<<"\n";
	while (sitr != send) {
		std::cerr<<std::setw(10)<<std::right<<sitr->first<<" "
				  <<std::setw(6)<<std::right<<sitr->second<<"\n";
		sitr++;
	}
}

}

}

#endif	/* TASKGENEGENEMODELREPORT_H */

