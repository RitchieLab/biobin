/* 
 * File:   taskbincollapse.h
 * Author: torstees
 *
 * Created on July 5, 2011, 12:12 PM
 * 
 * This structure represents a two stage model: 1) it builds the association between
 * groups and bins and 2) it dumps each individual's data to the binfile. 
 * 
 * During construction of the group/bin associations, a text report showing
 * the hierarchical relationship between a group and it's child group/regions/snps
 * can be written. 
 * 
 * For circular references (such as seen in go), the parent group is removed from
 * it's progeny's children-so, the SNPs found in the parent will not appear listed
 * under the child's SNPs. The only place I can see a problem is where Group A
 * fails to pass some threshold (Biofilter or Biobin) and contains B which also
 * contains A. If A contains any SNPs itself, they will not be expressed-neither
 * through A nor B. If A does pass the threshold and is collapsed (or used to
 * generate models in biofilter) everything works fine-since the SNPs in B's circular
 * reference of A are already accounted for in A itself....
 * 
 * Bin Data -- I'm keeping the bin output here in the collapse structure because
 * I don't want to have to carry that data around for extended periods. Also, 
 * the bin data within the individual is atomic in natural and can be used to 
 * form the larger entities-so this is plenty flexible enough. So, the knowledge
 * based bin information will only exist during the file production.
 * 
 */

#ifndef TASKBINCOLLAPSE_H
#define	TASKBINCOLLAPSE_H

#include <map>
#include "binapplication.h"
#include "biofilter/task.h"
#include "utility/locus.h"
#include <iostream>
namespace BioBin {
namespace Task {
	
class BinCollapse : public Biofilter::Task::Task {
public:
	BinCollapse();
	virtual ~BinCollapse();

	virtual void Init(Biofilter::Application* app);
	void ExecuteTask();
	void CollectSNPs(std::map<uint, 
			std::set<uint> >& snpCollection, 
			uint idx, 
			Knowledge::GroupManager& mgr, 
			Knowledge::RegionManager& regions);
	void WriteBinData();
	
	static int maxSnpCount;				///< This represents the threshold where we will drill down further in an attempt to find properly sized bins
	static bool VisualizeGroupTrees;
	static bool WriteKnowledgeBins;
protected:
	BinApplication *app;
	std::multimap<uint, uint> binContents;
	std::map<uint, uint> binIndex;
	void EvaluateGroup(uint depth, 
			std::map<uint, 
			std::set<uint> >& snpCollection, 
			std::set<uint>& visited, 
			uint idx, 
			Knowledge::GroupManager& gm, 
			Knowledge::RegionManager& regions);
	ofstream visualization;
	
	/**
	 * Associates a group with it's set of bin IDs. We can't store integers here, since 
	 * this is potentially related to more than one group manager
	 */
	std::vector<std::pair<Knowledge::Group*, Utility::IdCollection> > groupParticipant;	 
	Utility::IdCollection mergedBins;			///< Merged into one or more group
};
	
	
inline
BinCollapse::BinCollapse() : Task(3), app(NULL) {}

inline
BinCollapse::~BinCollapse() { }

inline
void BinCollapse::Init(Biofilter::Application* app) {
	this->app = (BinApplication*)app;
	
	
}

}
}

#endif	/* TASKBINCOLLAPSE_H */

