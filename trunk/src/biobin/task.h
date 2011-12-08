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

//#include "utility/types.h"
//#include "knowledge/regionmanager.h"
//#include "knowledge/snpdataset.h"
#include "application.h"

namespace BioBin {
namespace Task {

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

	virtual void Init(BioBin::Application* app);

	uint taskType;
	static bool detailedReport;
protected:
	std::string filename;
};

}

}

namespace std{

template <>
struct less<BioBin::Task::Task*>{
	bool operator()(const BioBin::Task::Task* x, const BioBin::Task::Task* y){
		return (x && y) ? (*x) < (*y) : y < x;
	}
};
}


using BioBin::Task::Task;

//typedef std::multimap<Task*, TaskPtrCmp> TaskList;
typedef std::multimap<uint, Task*> TaskList;
typedef std::pair<uint, Task*> TaskPair;

#endif	/* TASK_H */

