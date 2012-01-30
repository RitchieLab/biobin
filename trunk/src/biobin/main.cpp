#include "main.h"

#include "config.h"

#include <iostream>

#include "Configuration.h"
#include "knowledge/Configuration.h"

#include "task.h"
#include "taskfilegeneration.h"
#include "dataimporter.h"

#include <boost/program_options.hpp>

namespace po=boost::program_options;

using po::value;

/**
 * The VCF Tools require this, but I don't want to use their main...
 */
namespace VCF {
std::ofstream LOG;
}

namespace BioBin {

string Main::c_vcf_file;
string Main::c_knowledge_file;
string Main::c_genome_build = "37";

vector<string> Main::c_custom_groups;
set<uint> Main::c_source_ids;

/*
void Main::LoadConfiguration(const char *cfgFilename) {
	cfg.SetValue("REPORT_PREFIX", Utility::ExtractBaseFilename(cfgFilename));
	cfg.Parse(cfgFilename);
}
*/

void Main::InitRegionData() {
	vector<string> missingAliases;
	vector<string> aliasList;

	app.LoadRegionData(missingAliases, aliasList);
}

void Main::InitGroupData() {
	app.LoadGroupDataByName(c_custom_groups, c_source_ids);
}

void Main::RunCommands() {
	VCF::LOG.open("vcf-responses.log");

	app.Init(c_knowledge_file, true);
	if (c_genome_build != "") {
		app.LoadBuildConverter(c_genome_build);
	}

	//Tasks that run before SNPs load (not sure what those would be)
	RunTasks(0);

	/**
	 * Do the SNP oriented stuff here
	 */

	// TODO: check for existence of the file here!
	if(c_vcf_file.size() == 0){
		std::cerr<<"No SNP dataset available. Unable to continue.\n";
		exit(1);
	}

	LoadSNPs();

	RunTasks(1);

	InitRegionData();

	RunTasks(2);

	InitGroupData();
	app.InitBins();

	RunTasks(3);

}

void Main::RunTasks(int level){
	multimap<int, Task::Task*>::iterator itr = _task_list.lower_bound(level);
	multimap<int, Task::Task*>::iterator end = _task_list.upper_bound(level);

	while (itr != end) {
		(*itr).second->ExecuteTask();
		++itr;
	}
}

void Main::initTasks(){
	BioBin::Task::Task *t = new BioBin::Task::GenerateFiles(&app);
	_task_list.insert(std::make_pair(t->getType(), t));
}



void Main::LoadSNPs() {
	vector<string> lostSnps;
	app.InitVcfDataset(c_genome_build, lostSnps);
	std::cerr<<lostSnps.size()<<" SNPs were not able to be found in the variations database.\n";
}

}


int main(int argc, char *argv[]) {
	std::string cfgFilename;

	po::options_description cmd("General Options");
	cmd.add_options()
				("help,h","Display help message")
				("version,v","Display version")
				("sample-config,S", "Print a sample configuration to the screen");

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
	}catch(...){
		std::cout << "Error processing command line arguments\n";
		std::cout << cmd;
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
	}catch(...){
		std::cout << "Error processing configuration file\n";
		BioBin::Configuration::printConfig(std::cout);
		Knowledge::Configuration::printConfig(std::cout);
		return 2;
	}

	Knowledge::Configuration::parseOptions(vm);
	BioBin::Configuration::parseOptions(vm);

	BioBin::Main *app = new BioBin::Main();					///<The application object

	app->initTasks();

	/*if (!app->ParseCmdLine(argc, argv)) {
		delete app;
		exit(1);
	}*/
	//Performs any commands
	try {
		app->RunCommands();
	}
	catch (std::exception& e) {
		BioBin::BinApplication::errorExit = true;
		std::cerr<<"\nError: \t"<<e.what()<<".  Unable to continue.\n";
	}

	delete app;

	return 0;
}
