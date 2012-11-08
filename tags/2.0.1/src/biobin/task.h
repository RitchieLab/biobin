/* 
 * File:   task.h
 * Author: torstees
 * Tasks allow us to put specific tasks into a single object to avoid
 * overly complex meta objects (like the application)
 * Created on April 12, 2011, 4:51 PM
 */

#ifndef BIOBIN_TASK_TASK_H
#define	BIOBIN_TASK_TASK_H
#include <string>

#include <sys/types.h>

namespace BioBin {

class Application;

namespace Task {

class Task {

public:
	/**
	 *
	 * @param type Gene/SNP/Model
	 */
	Task(int taskType, BioBin::Application* app);

	virtual ~Task(){}
	virtual void ExecuteTask() = 0;

	//Required for STL container
	//bool operator<(const Task& other) const;

	virtual std::string GetFileSuffix() const { return ""; }
	virtual std::string GetFileExtension() const { return ""; }
	virtual std::string GetFileDescription() const { return "";}

	int getType() const {return _task_type;}

	//virtual void Init(BioBin::Application* app);

	static bool detailedReport;
protected:
	int _task_type;
	std::string _filename;
};

}

}
/*
namespace std{

template <>
struct less<BioBin::Task::Task*>{
	bool operator()(const BioBin::Task::Task* x, const BioBin::Task::Task* y){
		return (x && y) ? (*x) < (*y) : y < x;
	}
};
}
*/

//using BioBin::Task::Task;

//typedef std::multimap<Task*, TaskPtrCmp> TaskList;
//typedef std::multimap<uint, Task*> TaskList;
//typedef std::pair<uint, Task*> TaskPair;

#endif	/* TASK_H */

