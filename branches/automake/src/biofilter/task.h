/* 
 * File:   task.h
 * Author: torstees
 * Tasks allow us to put specific tasks into a single object to avoid
 * overly complex meta objects (like the application)
 * Created on April 12, 2011, 4:51 PM
 */

#ifndef TASK_H
#define	TASK_H
#include <map>

#include "utility/types.h"
#include "knowledge/regionmanager.h"
#include "knowledge/snpdataset.h"
#include "application.h"

namespace Biofilter {
namespace Task {
	class Task;
struct TaskPtrCmp
{
  bool operator()(const Task *t1, const Task *t2) const;
};

class Task;


class Task {
public:
	Task();
	/**
	 *
	 * @param type Gene/SNP/Model
	 */
	Task(uint taskType);

	virtual ~Task();

	virtual void ExecuteTask() = 0;

	//Required for STL container
	bool operator<(const Task& other) const;

	virtual std::string GetFileSuffix() { return ""; }
	virtual std::string GetFileExtension() { return ""; }
	virtual std::string GetFileDescription() { return "";}

	virtual void Init(Biofilter::Application* app);

	uint taskType;
	static bool detailedReport;
protected:
	Knowledge::RegionManager *regions;
	Knowledge::SnpDataset *snps;
	std::string filename;
};

inline
bool TaskPtrCmp::operator()(const Task *t1, const Task *t2) const {
	  return t1->taskType < t2->taskType;
 }

inline Task::Task() : regions(NULL), snps(NULL) { }

inline Task::Task(uint t) : taskType(t), regions(NULL), snps(NULL) { }


inline Task::~Task() {}
inline void Task::Init(Application* app) {
	regions			= app->GetRegions();
	snps				= app->GetDataset();
	filename			= app->AddReport(GetFileSuffix().c_str(), GetFileExtension().c_str(), GetFileDescription().c_str());
}

inline bool Task::operator<(const Task& other) const {
	return taskType<other.taskType;
}



}


}

//typedef std::multimap<Task*, TaskPtrCmp> TaskList;
typedef std::multimap<uint, Biofilter::Task::Task*> TaskList;
typedef std::pair<uint, Biofilter::Task::Task*> TaskPair;

#endif	/* TASK_H */

