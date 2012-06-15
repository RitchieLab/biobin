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

/*!
 * \brief The class that controls all top-level behavior.
 * This class, which should be used as a singleton (not enforced) controls all
 * of the top level behavior in BioBin.  It is a place that initializes all
 * of the data and runs the tasks that are given.
 *
 * Note that configuration processing is done before this class is instantiated
 * and is done in the main(int, char**) function.
 */
class Main {
public:

	/*!
	 * \brief Default constructor
	 * Creates an instance of the Main class.  Note: at this point, all
	 * configuration has been processed, so we can use the vcf file name here.
	 */
	Main() : app(c_knowledge_file, c_vcf_file){}
	~Main();

	/*!
	 * \brief Runs the commands in order.
	 * This is the main running function of the Main class.  Within this
	 * function, the tasks are run as follows:
	 * 0 - Run any setup tasks
	 * Load SNPs / Variants
	 * 1 - Run any tasks on only SNPS
	 * Load Region / gene data
	 * 2 - Run any tasks on Gene / SNP info only
	 * Load Pathway data
	 * Construct Bins
	 * 3 - Run any final tasks
	 */
	void RunCommands();

	/*!
	 * \brief Initializes tasks to be performed by the RunCommands function.
	 * This function sets up any tasks that need to be performed, and it must
	 * be called before RunComands.
	 *
	 * NOTE:  Should we make this private and call it from within RunCommands??
	 */
	void initTasks();

	//! The name of the VCF file containing info on variants
	static string c_vcf_file;
	//! The name of the SQLite LOKI database
	static string c_knowledge_file;
	//! The build of the genome that the VCF file is based on
	static string c_genome_build;

	//! A vector of custom groups to use
	static vector<string> c_custom_groups;
	//! A list of IDs of sources to use
	static set<uint> c_source_ids;

private:
	void LoadSNPs();
	void RunTasks(int level);
	void InitGroupData();
	void InitRegionData();
	
	BinApplication app;					///< The application that does all of the work

	multimap<int,Task::Task*> _task_list;

};

}

#endif
