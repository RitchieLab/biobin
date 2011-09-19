//
// C++ Implementation: appconfiguration
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "configuration.h"
#include <iomanip>
#include "knowledge/genegenemodel.h"
#include "knowledge/group.h"
#include "knowledge/region.h"
#include "knowledge/snpdataset.h"
#include "application.h"

#include "taskgenereport.h"
#include "taskmarkerinfo.h"
#include "tasksnpgenemap.h"
#include "tasksnpreport.h"
#include "taskgenegenemodelreport.h"
#include "tasksnpsnpmodelarchive.h"

#ifdef LOCAL_RELEASE
#define BIODB "/projects/ritchie/biofilter/knowledge.bio"
#else
#define BIODB "knowledge.bio"
#endif
using namespace std;

namespace Biofilter {

Configuration::Configuration() {}


Configuration::~Configuration() {
	TaskList::iterator itr = tasks.begin();
	TaskList::iterator end = tasks.end();

	while (itr != end) 
		delete itr++->second;

}

void Configuration::Init() {
	InitKey("SETTINGS_DB",				Utility::ENV("SETTINGS_DB", BIODB).c_str(),  "BioFilter data");
	InitKey("MAX_GENE_COUNT",			"30",					"Max number of genes before we ignore the group");
	InitKey("RS_SOURCE",					"",					"The source file for the RS numbers in your dataset\n#This is used as an alternative to a MAP_SOURCE file");
	InitKey("MAP_SOURCE",				"",					"User specified chrom-rsid-gen_dist-bp_location (plink 4 column format)\n#Chromosomes can be 1-22,X,Y,XY,MT\n#Must be accompanied with MAP_GENOME_BUILD");
	InitKey("ADD_GROUP",					"",					"Add user defined group-there can be more than one ADD_GROUP line in a single configuration");
	//	InitKey("MAP_GENOME_BUILD",		"",					"Specify the genomic build associated with MAP_SOURCE.\n#This isn't used with RS_SOURCE");
	InitKey("INCLUDE_GROUPS",			"",					"List the various groups (by group ID) separated by spaces");
	InitKey("INCLUDE_GROUP_FILE",		"",					"File containing group IDs to be the groups to be searched");
	InitKey("INCLUDE_GROUP_NAMES",	"",					"List various groups (by name) separated by spaces. The name most be spelled EXACTLY as it is in the database.");
	InitKey("INCLUDE_GROUP_NAME_FILE", "",					"File containing group Names to be part of the search.");
	//InitKey("MODEL_FILENAME", 		"NONE",				"Set the filename for the output model list (none writes to std-out)");
//	InitKey("MODEL_BUFFER_INIT",		"1000000",			"Set the initial size of the model buffer. ");
//	InitKey("MODEL_BUFFER_MAX",		"10000000",			"Set the upper limit to the buffer. Bigger -> faster, but must remain within\n# the limits of the hardware or could cause the application\n# to fail or become so slow that it will never complete.");
	InitKey("POPULATION",				"NO-LD",				"Set the population ID to match the population your data is drawn from so that\n# LD patterns can be used to expand the gene boundaries.");
	InitKey("GENE_BOUNDARY_EXTENSION", "0",			"How many base pair locations up and down stream do we expand gene boundaries (Only used if POPULATION is NO-LD)");
	//	InitKey("DISEASE_DEPENDENT",		"",					"Add one or more files containing disease dependent genes ");
//	InitKey("PREFERRED_ALIAS",			"",					"User can specify aliases for genes (the alias must be present in the database");
	InitKey("REPORT_PREFIX", 			"",		 			"	Prefix used for all reports");
//	InitKey("LOAD_ALL_ALIASES",		"NO",					"Loads all aliases and generates a text report containing their associations");
//	InitKey("HTML_REPORTS",				"NO",					"Write reports in html format (not all reports support HTML formatting");
	InitKey("IMPLICATION_IDX_DUPLICATE_WEIGHT", "0.0",	"Weight applied to implication index for disease dependent groups are associated with both genes");
	InitKey("BINARY_MODEL_ARCHIVE",	"YES",				"Indicates whether to use a binary format instead of a text version");
	InitKey("DISEASE_DEPENDENT_LEVEL", "ALL_MODELS",	"ALL_MODELS, GROUP_LEVEL, DD_ONLY  These are used to determine selectivity of gene/gene models based on disease dependent relationships.");
//	InitKey("COLLAPSE_ASSOCIATION_REPORT", "NO",			"When true, the associations reported cease where groups would generate models (if possible)");
//	InitKey("ASSOCIATION_REPORT",		"NO",					"Produces association report");
//	InitKey("ASSOCIATION_GRAPH",		"NO",					"Produces association graph input files");
//	InitKey("CLEANUP_RSIDS",			"ON",					"Attempt to identify merged and expired RS IDs in the SNP Source and modify them accordingly");
	InitKey("GENOMIC_BUILD",			"37",					"Determine what build any map files are based on.");
	InitKey("DETAILED_REPORTS",		"OFF",				"Activates extra data in the reports. See manual for details.");
	InitKey("MARKER_INFO_REPORT",		"OFF",				"Produces a marker-info report once the SNP data is loaded.");
	InitKey("SNP_REPORT",				"OFF",				"Produces comma separated list of genes each SNP is found in.");
	InitKey("SNP_GENE_MAP",				"OFF",				"Produces report indicating SNP/gene relationship including relational details (Interior/etc)");
	InitKey("GENE_COVERAGE",			"",					"Filename (or ALL) containing list of gene aliases of to be used in various gene reports.");
	InitKey("GENE_REPORT",				"OFF",				"Details the contents of the regions in use.");
	InitKey("COVERAGE_RS",				"",					"Platforms for coverage reports using RS IDs instead of MAP files");
	InitKey("COVERAGE_MAP",				"",					"Platforms for coverage reports using MAP files instead of RS IDs");

	InitKey("MINIMUM_IMPLICATION_INDEX", "-1",			"The minimum implication index to show in gene/gene or SNP/SNP models");
	InitKey("MAX_SNP_MODEL_COUNT",			"-1",					"The max number of SNP/SNP models to be generated");
	InitKey("EXPORT_SNP_MODELS",		"NO",					"Exports snp-snp models immediately after gene-gene model production");
	InitKey("EXPORT_GENE_MODELS",		"NO",					"Generate gene-gene models");

	/**
	 */
}
void Configuration::PrintSet(const char *key, vector<string>& settings, ostream& os) {
	vector<string>::iterator itr = settings.begin();
	vector<string>::iterator end = settings.end();

	os<<setw(35)<<right<<key<<" : ";
	int count = 0;
	while (itr != end ) {
		if (count++>0)
			os<<",";
		os<<*itr++;
	}
	os<<"\n";
}

void Configuration::ReportConfiguration(std::ostream& os) {
	map<string, vector<string> >::iterator itr= strings.begin();
	map<string, vector<string> >::iterator end = strings.end();

	os<<"-------------------- Configuration Parameters ----------\n";
	while (itr != end) {
		PrintSet(itr->first.c_str(), itr->second, os);
		itr++;
	}	
}

void Configuration::WriteConfiguration(std::ostream& os) {

}

void Configuration::LoadFileContents(const char *key, Utility::IdCollection& fileContents) {

	Utility::StringArray filenames;
	GetLines(key, filenames);

	Utility::StringArray::iterator itr = filenames.begin();
	Utility::StringArray::iterator end = filenames.end();

	while (itr != end) {
		Utility::IdCollection ids;
		std::string contents = Utility::LoadContents(itr->c_str());
		ids = Utility::ToSet<uint>(contents.c_str(), "\n");
		fileContents.insert(ids.begin(), ids.end());
	}

}

void Configuration::LoadFileContents(const char *key, Utility::StringArray& fileContents) {

	Utility::StringArray filenames;
	GetLines(key, filenames);

	if (filenames.size() == 1 && filenames[0] == "ALL") {
		fileContents.push_back("ALL");
		return;
	}

	Utility::StringArray::iterator itr = filenames.begin();
	Utility::StringArray::iterator end = filenames.end();

	while (itr != end) {
		std::string contents = Utility::LoadContents(itr++->c_str());
		Utility::StringArray lines = Utility::Split(contents.c_str(), "\n");
		fileContents.insert(fileContents.end(), lines.begin(), lines.end());
	}
	
}

int Configuration::RunTasks(uint taskType) {
	TaskList::iterator itr = tasks.lower_bound(taskType);
	TaskList::iterator end = tasks.upper_bound(taskType);

	int count = 0;
	while (itr != end) {
		(itr++)->second->ExecuteTask();
		count++;
	}
	return count;
}

int Configuration::CountTasks(uint taskType) {
	return tasks.count(taskType);
}

void Configuration::AddTask(const char *key, Task::Task* item) {
	if (GetBoolean(key)) 
		tasks.insert(TaskPair(item->taskType, item));
	else
		delete item;
}

void Configuration::ExecuteConfiguration(Application* app) {
	//Write all of our settings to the relevant variables in memory
	Knowledge::Region::DuplicateDD_Weight				= GetDouble("IMPLICATION_IDX_DUPLICATE_WEIGHT");
	Knowledge::RegionManager::modelGenerationType	= Knowledge::ModelGenerationMode::ConvertType(GetString("DISEASE_DEPENDENT_LEVEL").c_str());
	Knowledge::BinaryArchive								= GetBoolean("BINARY_MODEL_ARCHIVE");
	Task::Task::detailedReport								= GetBoolean("DETAILED_REPORTS");
	Knowledge::GeneGeneModelArchive::minImplicationIndex		= GetInteger("MINIMUM_IMPLICATION_INDEX");
	Knowledge::GeneGeneModelArchive::maxModelCount	= GetInteger("MAX_SNP_MODEL_COUNT");


	//Knowledge::KbGroup::CollapseAssociationReport = GetBoolean("COLLAPSE_ASSOCIATION_REPORT");
	app->SetGeneExtension(GetInteger("GENE_BOUNDARY_EXTENSION"));
	app->SetReportPrefix(GetString("REPORT_PREFIX").c_str());
	app->UseHtmlReports(GetBoolean("HTML_REPORTS"));
	
	


	//Build out the task list
	AddTask("GENE_REPORT", new Task::GeneReport());
	AddTask("MARKER_INFO_REPORT", new Task::MarkerInfo());
	AddTask("SNP_GENE_MAP", new Task::SnpReport());
	AddTask("SNP_REPORT", new Task::SnpGeneMap());
	AddTask("EXPORT_GENE_MODELS", new Task::GeneGeneModelReport());
	AddTask("EXPORT_SNP_MODELS", new Task::SnpSnpModelArchive());
	
	TaskList::iterator itr = tasks.begin();
	TaskList::iterator end = tasks.end();

	while (itr != end)
		itr++->second->Init(app);

	
}

}
