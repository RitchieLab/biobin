#ifndef MAIN_H
#define MAIN_H

//#include "configuration.h"

#include <string>
#include <vector>
#include <set>
#include <sys/types.h>

#include "binapplication.h"

using std::set;
using std::vector;
using std::string;

namespace BioBin {

namespace Task{
class Task;
}

class Main {
public:
	Main() : app(c_vcf_file){}
	~Main(){}
	/**
	 * @brief Pass the arguments to the application object
	 */
	//bool ParseCmdLine(int argc, char **argv);
	//int ParseCmd(int curr, int argc, char **argv);
	//void PrintHelp();						///< Display usage details
	//void PrintBanner();					///< Display details about the software
	//void LoadConfiguration(const char *cfgFilename);
	void RunCommands();
	//int SetConfigValue(int nextCmd, int argc, const char *var, const char *val, const char *err);
	void InitGroupData();
	void InitRegionData();

	void initTasks();

	void RunTasks(int level);

	static string c_vcf_file;
	static string c_knowledge_file;
	static string c_genome_build;

	static vector<string> c_custom_groups;
	static set<uint> c_source_ids;

protected:
	void LoadSNPs();

	
	//Configuration cfg;					///< Configuration settings
	BinApplication app;					///< The application that does all of the work
	//bool silentRun;						///< Used to silence the banner banter

	multimap<int,Task::Task*> _task_list;
};

}

#endif
