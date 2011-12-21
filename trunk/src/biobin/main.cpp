#include "main.h"

#include "config.h"

#include <iostream>
#include "utility/exception.h"

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

	DataImporter vcfimporter(c_vcf_file);
	LoadSNPs(vcfimporter);

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

/*

bool Main::ParseCmdLine(int argc, char **argv) {

	//Test the DB connection
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
	if (argc < 2) {
		PrintHelp();
		return false;
	}
	int i=1;
	cfg.Init();
	if (argv[1][0] != '-')
		LoadConfiguration(argv[i++]);
	//Work out any other cmd line arguments
	for (; i<argc && i>0;) {
		i=ParseCmd(i, argc, argv);
	}
	cfg.ExecuteConfiguration(&app);
	app.SetReportPrefix(cfg.GetLine("REPORT_PREFIX").c_str());

	if (action == BiofilterAction::ParseError) {
		return false;
	}
	if (action == BiofilterAction::PrintSampleConfig) {
		PrintBanner();
		std::cout<<"#BioBin configuration file\n";
		std::cout<<"#\n#Users can change these parameters to meet their needs.\n";
		std::cout<<"#Please see the manual for more information about the different parameters and their options.\n";
		cfg.Write(std::cout);
		return false;
	}

	if (!silentRun)
		cfg.ReportConfiguration(std::cerr);

	return true;
}
*/

void Main::LoadSNPs(DataImporter& vcf) {
	vector<string> lostSnps;
	app.InitVcfDataset(c_genome_build, lostSnps, vcf);
	std::cerr<<lostSnps.size()<<" SNPs were not able to be found in the variations database.\n";
}

/*
int Main::SetConfigValue(int nextCmd, int argc, const char *var, const char *val, const char *err) {
	if (nextCmd < argc) {
		cfg.SetValue(var, val);
	} else {
		action = BiofilterAction::ParseError;
		std::cerr<<err<<"\n";
		return -1;
	}
	return nextCmd + 1;
}

int Main::ParseCmd(int curr, int argc, char **argv) {
	int nextCmd = curr+1;
	if (strcmp(argv[curr], "-S")==0 || strcmp(argv[curr], "--sample-config")==0) {
		action = BiofilterAction::PrintSampleConfig;
		return nextCmd;
	}
	if (strcmp(argv[curr], "--DB")==0){
		return SetConfigValue(nextCmd, argc, "SETTINGS_DB", argv[nextCmd], "--DB must be followed by a database filename");
	}
	if (strcmp(argv[curr], "-b")==0 || strcmp(argv[curr], "--binary")==0){
		return SetConfigValue(nextCmd, argc, "BINARY_MODEL_ARCHIVE", argv[nextCmd], "--binary must be followed by Yes/No");
	}
	if (strcmp(argv[curr], "-D")==0) {
		cfg.SetValue("DETAILED_REPORTS", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "-B")==0 || strcmp(argv[curr], "--build")==0) {
		if (nextCmd < argc) {
			cfg.SetValue("GENOMIC_BUILD", argv[nextCmd++]);
			return nextCmd;
		} else {
			action = BiofilterAction::ParseError;
			std::cerr<<"--build must be followed by an appropriate build number (35, 36, etc.)\n";
			return -1;
		}
	}

	if (strcmp(argv[curr], "-k")==0 || strcmp(argv[curr], "--knowledge-threshold")==0){
		return SetConfigValue(nextCmd, argc, "BIN_COLLAPSE_THRESHOLD", argv[nextCmd], "--knowledge-threshold must be followed by the max number of SNPs allowed in knowledge based bins.");
	}
	if (strcmp(argv[curr], "-t")==0 || strcmp(argv[curr], "--maf-threshold")==0){
		return SetConfigValue(nextCmd, argc, "MAF_CUTOFF", argv[nextCmd], "--MAF_CUTOFF must be followed by maf threshold value");
	}
	if (strcmp(argv[curr], "--PREFIX")==0){
		return SetConfigValue(nextCmd, argc, "REPORT_PREFIX", argv[nextCmd], "--PREFIX must be followed by prefix to be prepended to the generated filenames");
	}
	if (strcmp(argv[curr], "-p")==0 || strcmp(argv[curr], "--set-population")==0){
		return SetConfigValue(nextCmd, argc, "POPULATION", argv[nextCmd], "--set-population must be followed by name population you wish to use");
	}
	if (strcmp(argv[curr], "--gene-boundary")==0){
		return SetConfigValue(nextCmd, argc, "GENE_BOUNDARY_EXTENSION", argv[nextCmd], "--gene-boundary must be followed by an integer describing the number of bases");
	}
	if (strcmp(argv[curr], "-V")==0 || strcmp(argv[curr], "--vcf-file")==0){
		return SetConfigValue(nextCmd, argc, "VCF_FILE", argv[nextCmd], "--vcf-file must be followed by the name of a vcf file.");
	}
	if (strcmp(argv[curr], "-z")==0 || strcmp(argv[curr], "--gzipped vcf")==0){
		return SetConfigValue(nextCmd, argc, "COMPRESSED_VCF", argv[nextCmd], "--gzipped vcf must be followed by YES/NO.");
	}
	action = BiofilterAction::ParseError;
	std::cerr<<"Unrecognized parameter: "<<argv[curr]<<"\n";
	return -1;
}
*/

}


int main(int argc, char *argv[]) {
	std::string cfgFilename;

	BioBin::Main *app = new BioBin::Main();					///<The application object

	po::options_description cmd("Biobin Options");
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

	app->initTasks();

	/*if (!app->ParseCmdLine(argc, argv)) {
		delete app;
		exit(1);
	}*/
	//Performs any commands
	try {
		app->RunCommands();
	}
	catch (Utility::Exception::General& e) {
		BioBin::BinApplication::errorExit = true;
		std::cerr<<"\nError: \t"<<e.GetErrorMessage()<<" Unable to continue.\n";
	}

	delete app;

	return 0;
}
