/* 
 * File:   task.cpp
 * Author: torstees
 * 
 * Created on April 12, 2011, 4:51 PM
 */

#include "task.h"

namespace BioBin {
namespace Task {

//Task::Task() : regions(NULL), snps(NULL) { }

Task::Task(uint t) : taskType(t) { }

Task::~Task() {}

void Task::Init(Application* app) {
	filename = app->AddReport(GetFileSuffix().c_str(), GetFileExtension().c_str(), GetFileDescription().c_str());
}

bool Task::operator<(const Task& other) const {
	return taskType<other.taskType;
}

bool Task::detailedReport = true;
}
}

