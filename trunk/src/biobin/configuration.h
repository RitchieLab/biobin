//
// C++ Interface: appconfiguration
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOBINAPPCONFIGURATION_H
#define BIOBINAPPCONFIGURATION_H

#include "biofilter/application.h"
#include "utility/lineparser.h"
#include "biofilter/task.h"
#include "binapplication.h"
namespace BioBin {


/**
@Brief Reads/writes Configuration and initializes application settings accordingly

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Configuration : public Utility::FileToMapExplicit {
public:
    Configuration();
    ~Configuration();

	/**
	 * @Brief Write the configuration to stream
	 */
	void WriteConfiguration(std::ostream& os);

	/**
	 * @Brief Initialize configuration with defaults
	 */
	void Init();
	 
	void AddTask(const char *key, Biofilter::Task::Task* item);

	int RunTasks(uint taskType);

	int CountTasks(uint taskType);
	
	/**
	 * @Brief Updates application data according to values from the configuration
	 */
	void ExecuteConfiguration(BinApplication* app);

	/**
	 * @Brief Produce a report header listing key settings from the configuration
	 */
	void ReportConfiguration(std::ostream& os);

	void PrintSet(const char *key, std::vector<std::string>& settings, std::ostream& os);

	void LoadFileContents(const char *key, Utility::StringArray& fileContents);

	void LoadFileContents(const char *key, Utility::IdCollection& fileContents);

	TaskList tasks;

};



}

#endif
