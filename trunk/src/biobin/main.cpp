#include "main.h"

#include "config.h"

#include <iostream>

#include "Configuration.h"
#include "knowledge/Configuration.h"
#include "binapplication.h"

#include "task.h"
#include "taskfilegeneration.h"

// Use the boost filesystem library to work with OS-independent paths
#include <boost/filesystem.hpp>

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

string Main::c_vcf_file="";
string Main::c_knowledge_file="";
string Main::c_genome_build = "37";

vector<string> Main::c_custom_groups;

Main::~Main(){
	multimap<int, Task::Task*>::iterator m_itr = _task_list.begin();
	while(m_itr != _task_list.end()){
		delete (*m_itr).second;
		++m_itr;
	}
	_task_list.clear();
}

void Main::initTasks(){
	BioBin::Task::Task *t = new BioBin::Task::GenerateFiles(&app);
	_task_list.insert(std::make_pair(t->getType(), t));
}

void Main::RunCommands() {

	//VCF::LOG.open("vcf-responses.log");

	//app.Init(c_knowledge_file, true);

	//Tasks that run before SNPs load (not sure what those would be)
	RunTasks(0);

	LoadSNPs();

	RunTasks(1);

	InitRegionData();

	RunTasks(2);

	InitGroupData();
	app.InitBins();

	RunTasks(3);

}

void Main::InitRegionData() {
	vector<string> missingAliases;
	vector<string> aliasList;

	app.LoadRegionData(missingAliases, aliasList);
}

void Main::InitGroupData() {
	app.LoadGroupDataByName(c_custom_groups);
}

void Main::RunTasks(int level){
	multimap<int, Task::Task*>::iterator itr = _task_list.lower_bound(level);
	multimap<int, Task::Task*>::iterator end = _task_list.upper_bound(level);

	while (itr != end) {
		(*itr).second->ExecuteTask();
		++itr;
	}
}





void Main::LoadSNPs() {
	vector<string> lostSnps;
	app.InitVcfDataset(c_genome_build, lostSnps);
	if (lostSnps.size() > 0){
		std::cerr << "WARNING: Could not translate " << lostSnps.size() << " variants.\n";
	}
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
	}catch(...){
		std::cout << "Error processing configuration file\n";
		BioBin::Configuration::printConfig(std::cout);
		Knowledge::Configuration::printConfig(std::cout);
		return 2;
	}

	try{
		Knowledge::Configuration::parseOptions(vm);
		BioBin::Configuration::parseOptions(vm);
	}catch(...){
		std::cerr<<"\nError Parsing Configuration File\n";
		return 3;
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
	}

	BioBin::Main *app = new BioBin::Main();					///<The application object

	if(BioBin::BinApplication::s_run_normal){
		app->initTasks();

		try {
			app->RunCommands();
		}
		catch (std::exception& e) {
			BioBin::BinApplication::errorExit = true;
			std::cerr<<"\nError: \t"<<e.what()<<".  Unable to continue.\n";
		}
	}

	delete app;

	return 0;
}
