/* 
 * File:   taskgenereport.h
 * Author: torstees
 *
 * Created on April 12, 2011, 5:04 PM
 */
#ifndef TASK_MARKER_INFO_H
#define	TASK_MARKER_INFO_H

#include "task.h"
#include <string>


namespace Biofilter {
namespace Task {

class MarkerInfo : public Task {
public:
	MarkerInfo();

	MarkerInfo(const MarkerInfo& orig);

	virtual ~MarkerInfo();

	virtual std::string GetFileSuffix() { return "marker-info"; }
	virtual std::string GetFileExtension() { return "map"; }
	virtual std::string GetFileDescription() { return "Marker Info Report";}
	void ExecuteTask();

};

inline
MarkerInfo::MarkerInfo() : Task(1) {}

inline
MarkerInfo::~MarkerInfo() { }

inline
void MarkerInfo::ExecuteTask() {
	snps->WriteMarkerInfo(filename.c_str());
}

}
}

#endif	/* TASKGENEREPORT_H */
