#include "main.h"

//#include "config.h"

#include <iostream>
#include <cstring>
#include <cstdlib>

#include "Configuration.h"
#include "knowledge/Configuration.h"
#include "binmanager.h"

// Use the boost filesystem library to work with OS-independent paths
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// We'll need this to set the GSL error handler
#include <gsl/gsl_errno.h>
#ifdef HAVE_EXECINFO
#include <execinfo.h>
#endif
#ifdef HAVE_GCC_ABI_DEMANGLE
#include <cxxabi.h>
#endif


namespace po=boost::program_options;

using po::value;
using std::string;
using std::vector;
using std::multimap;
using std::ifstream;
using std::cerr;


namespace BioBin {

string Main::c_vcf_file="";
string Main::c_knowledge_file="";
string Main::c_genome_build = "37";
string Main::OutputDelimiter = ",";

bool Main::WriteBinData = true;
bool Main::WriteLociData = true;


vector<string> Main::c_custom_groups;

Main::~Main(){}

void Main::RunCommands() {

	app.InitVcfDataset(c_genome_build);

	vector<string> missingAliases;
	vector<string> aliasList;

	app.LoadRegionData(missingAliases, aliasList);

	// only do this if we're expanding by role, please!
	if(BinManager::ExpandByExons){
		app.loadRoles();
	}

	// only do this if we're binning by pathway, please!
	if(BinManager::UsePathways){
		app.LoadGroupDataByName(c_custom_groups);
	}

	app.InitBins();

	if (WriteLociData){
		std::string filename = app.reportPrefix + "-locus.csv";
		app.writeLoci(filename,OutputDelimiter);
	}

}

void Main::gsl_tracer(const char* reason, const char* filename, int line, int gsl_error){
	cerr<<"WARNING: GSL Error " << gsl_error << " detected:" << std::endl;
	cerr << reason << " in " << filename << ", line " << line << std::endl;
#ifdef HAVE_EXECINFO
	const unsigned int buf_size=50;
	void* buffer[buf_size];
	int n_buf = backtrace(buffer, buf_size);
	cerr << "Backtrace (maximum depth of " << buf_size << "):" << std::endl;
	char** ptr = backtrace_symbols(buffer, n_buf);

	for (int idx = 0; idx < n_buf; idx++) {
		cerr << idx << ":";
#ifdef HAVE_GCC_ABI_DEMANGLE
		// Here, we're able to demangle the names, so let's do that!
		// first, find the position of the first open paren
		char* tok = strchr(ptr[idx], '(');
		int paren_pos = tok - ptr[idx];

		// Now, find the position of the first "+" to occur after the open paren
		tok = strchr(ptr[idx] + (paren_pos + 1), '+');
		char* paren_tok = strchr(ptr[idx] + (paren_pos + 1), ')');

		if (tok && paren_tok > tok) {

			// convert the open paren to a NUL
			ptr[idx][paren_pos] = '\0';

			int plus_pos = tok - ptr[idx];

			// convert that plus sign to a NUL
			ptr[idx][plus_pos] = '\0';

			// great! now, we can get that string and demangle it
			char* demangled_name = NULL;
			int status;
			size_t demangled_buflen;

			demangled_name = abi::__cxa_demangle(ptr[idx] + (paren_pos + 1),
					demangled_name, &demangled_buflen, &status);

			if (status == 0) {
				// print up to (not including) the first paren, the demangled name, then everything after the '+'
				cerr << ptr[idx] << "(" << demangled_name << "+"
						<< ptr[idx] + (plus_pos+1);
				free(demangled_name);
			} else {
				cerr << ptr[idx] << "(" << ptr[idx] + (paren_pos + 1) << "+"
						<< ptr[idx] + (plus_pos + 1);
			}

		} else {
			cerr << ptr[idx];
		}
#else
		cerr << ptr[idx];
#endif
		cerr << std::endl;

	}

	free(ptr);
#else
	cerr << "No backtrace available" << std::endl;
#endif
}

}


int main(int argc, char *argv[]) {
	std::string cfgFilename;

	po::options_description cmd("General Options");
	cmd.add_options()
				("help,h","Display help message")
				("version,v","Display version")
				("sample-config,S", "Print a sample configuration to the screen")
				("no-gsl-bt", "Abort on GSL Errors; do not try to print backtrace");

	Knowledge::Configuration::addCmdLine(BioBin::Configuration::addCmdLine(cmd));

	po::options_description hidden("Hidden Biobin Options");
	hidden.add_options()
					("config-file",value<vector<string> >(),"Name of the configuration file");

	po::options_description cmd_options;
	cmd_options.add(cmd).add(hidden);

	po::options_description config_options;
	Knowledge::Configuration::addConfigFile(BioBin::Configuration::addConfigFile(config_options));

	po::positional_options_description pos;
	pos.add("config-file",-1);

	po::variables_map vm;
	try{
		store(po::command_line_parser(argc,argv).options(cmd_options).positional(pos).run(), vm);
		notify(vm);
	}catch(std::exception& e){
		std::cerr << "Error processing command line arguments, please see the --help option for more details\n";
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	// Check here for help, version or sample printing
	if (vm.count("help")){
		std::cout << cmd;
		return 1;
	}

	if (vm.count("sample-config")){
		BioBin::Configuration::printConfig(std::cout);
		Knowledge::Configuration::printConfig(std::cout);
		return 1;
	}

	if (vm.count("version")){
		std::cout << PACKAGE_STRING << "\n";
		std::cout << "(c) Ritchie Lab, 2012\n";
		std::cout << "To report bugs, please email " << PACKAGE_BUGREPORT << "\n";
		BioBin::BinApplication::s_run_normal = false;
	}

	try{
		// load the info given in the configuration file
		if(vm.count("config-file")){
			vector<string> conf_files = vm["config-file"].as<vector<string> >();
			vector<string>::const_iterator itr = conf_files.begin();
			vector<string>::const_iterator end = conf_files.end();
			while(itr != end){
				ifstream ifs((*itr).c_str());
				if (!ifs)
				{
					std::cerr << "WARNING: can not open config file: " << (*itr) << "\n";
				}
				else
				{
					store(parse_config_file(ifs, config_options), vm);
					notify(vm);
				}
				++itr;
			}
		}
	}//catch(...){
	 catch(po::error& e){
		BioBin::Configuration::printConfig(std::cout);
		Knowledge::Configuration::printConfig(std::cout);
		std::cout << "\n#### Error processing configuration file ####\n";
		std::cout << e.what() << std::endl;
		return 2;
	}catch(...){
                BioBin::Configuration::printConfig(std::cout);
                Knowledge::Configuration::printConfig(std::cout);
		std::cout << "\n#### Error processing configuration file ####\n";
        }
	

	try{
		Knowledge::Configuration::parseOptions(vm);
		BioBin::Configuration::parseOptions(vm);
	}catch(...){
		std::cerr<<"\nError Parsing Configuration File\n";
		return 3;
	}

	gsl_error_handler_t *old_handler = NULL;
	if(!vm.count("no-gsl-bt")){
		// first things first, let's set that GSL handler!
		old_handler = gsl_set_error_handler(&BioBin::Main::gsl_tracer);

	}


	if(BioBin::BinApplication::s_run_normal){

		// TODO: check for existence of the file here!
		if(BioBin::Main::c_vcf_file.size() == 0){
			std::cerr<<"ERROR: No VCF file given.  You must supply a vcf file.\n";
			exit(1);
		}

		boost::filesystem::path vcf_path = boost::filesystem::path(BioBin::Main::c_vcf_file);
		if (!boost::filesystem::is_regular_file(vcf_path)) {
			std::cerr<<"ERROR: Could not find VCF file at " << vcf_path << "\n";
			exit(1);
		}

		boost::filesystem::path report_path=boost::filesystem::absolute(BioBin::BinApplication::reportPrefix).parent_path();

		if( !boost::filesystem::is_directory(report_path) ){
			std::cerr<<"WARNING: report-prefix path does not exist, attempting to create" << std::endl;
			boost::system::error_code ec;
			if(!boost::filesystem::create_directories(report_path, ec)){
				std::cerr << "ERROR: could not create directory given in report-prefix" << std::endl;
				exit(1);
			}
		}
	}

	BioBin::Main *app = new BioBin::Main();					///<The application object

	if(BioBin::BinApplication::s_run_normal){

		try {
			app->RunCommands();
		}
		catch (std::exception& e) {
			BioBin::BinApplication::errorExit = true;
			std::cerr<<"\nError: \t"<<e.what()<<".  Unable to continue.\n";
		}
	}

	delete app;

	if(!vm.count("no-gsl-bt")){
		// unset the GSL handler here, please!
		gsl_set_error_handler(old_handler);
	}

	return 0;
}
