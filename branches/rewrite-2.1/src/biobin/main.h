#ifndef MAIN_H
#define MAIN_H

//#include "configuration.h"

#include <string>
#include <vector>
#include <set>
#include <sys/types.h>

#include "binapplication.h"

namespace BioBin {


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

	// does the actual biobin work (only called if configuration is A-OK
	void RunCommands();

	//! The name of the VCF file containing info on variants
	static std::string c_vcf_file;
	//! The name of the SQLite LOKI database
	static std::string c_knowledge_file;
	//! The build of the genome that the VCF file is based on
	static std::string c_genome_build;

	static std::string OutputDelimiter;

	static bool WriteLociData;
	static bool WriteBinData;

	//! A vector of custom groups to use
	static std::vector<std::string> c_custom_groups;


private:
	
	BinApplication app;					///< The application that does all of the work

};

}

#endif
