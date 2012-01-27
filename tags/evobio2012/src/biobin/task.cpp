/* 
 * File:   task.cpp
 * Author: torstees
 * 
 * Created on April 12, 2011, 4:51 PM
 */

#include "task.h"

#include "application.h"

namespace BioBin {
namespace Task {

bool Task::detailedReport = true;

Task::Task(int t, Application* app) : _task_type(t) {
	app->AddReport(GetFileSuffix(), GetFileExtension(), GetFileDescription());
}

}
}

