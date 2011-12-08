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
//#include "knowledge/genegenemodel.h"
//#include "knowledge/group.h"
//#include "knowledge/region.h"
//#include "knowledge/snpdataset.h"

#include "taskfilegeneration.h"
//#include "biofilter/taskgenereport.h"
//#include "biofilter/taskmarkerinfo.h"
//#include "biofilter/tasksnpgenemap.h"
//#include "biofilter/tasksnpreport.h"
//#include "biofilter/taskgenegenemodelreport.h"
//#include "biofilter/tasksnpsnpmodelarchive.h"
#include <boost/algorithm/string.hpp>
#include "binmanager.h"
#include "dataimporter.h"
//#include "taskbincollapse.h"

#ifdef LOCAL_RELEASE
#define BIODB "/projects/ritchie/knowledge.bio"
#else
#define BIODB "knowledge.bio"
#endif
using namespace std;

namespace BioBin {

Configuration::Configuration() {}


Configuration::~Configuration() {
	/*
	TaskList::iterator itr = tasks.begin();
	TaskList::iterator end = tasks.end();

	while (itr != end) 
		delete itr++->second;
	 */
}

void Configuration::Init() {
	InitKey("SETTINGS_DB",				Utility::ENV("SETTINGS_DB", BIODB).c_str(),  "BioFilter data");
	InitKey("VCF_FILE",				"",					"List of vcf files associated with the data. Currently, biobin assumes a SNP occurs in only one file.");
	InitKey("COMPRESSED_VCF",		"NO",					"YES/NO gzipped VCF files.");
	InitKey("MAF_CUTOFF",				"0.05",			"Threshold associated with calling rare-variants.");
	InitKey("INCLUDE_GROUPS",			"",				"List the various groups (by group ID) separated by spaces");
	InitKey("INCLUDE_GROUP_FILE",		"",				"File containing group IDs to be the groups to be searched");
	InitKey("INCLUDE_GROUP_NAMES",	"",				"List various groups (by name) separated by spaces. The name most be spelled EXACTLY as it is in the database.");
	//InitKey("INCLUDE_GROUP_NAME_FILE", "",					"File containing group Names to be part of the search.");
	InitKey("POPULATION",				"NO-LD",		   "Set the population ID to match the population your data is drawn from so that\n# LD patterns can be used to expand the gene boundaries.");
	InitKey("GENE_BOUNDARY_EXTENSION", "0",			"How many base pair locations up and down stream do we expand gene boundaries (Only used if POPULATION is NO-LD)");
	//	InitKey("DISEASE_DEPENDENT",		"",					"Add one or more files containing disease dependent genes ");
	//InitKey("REPORT_PREFIX", 			"",		 			"	Prefix used for all reports");
	InitKey("GENOMIC_BUILD",			"37",			   "Determine what build any map files are based on.");
	InitKey("WRITE_BIN_DATA",        "YES",         "Writes the bin counts to the file");
	InitKey("WRITE_GENOTYPE_DATA",   "YES",         "Writes the genotype counts to the file");
	InitSingletary("OUTPUT_DELIMITER",      "'	'","The string to be used to delimit fields in the data output (defaut is a single space).");
	InitKey("PHENOTYPE_FILENAME",    "",            "Phenotype file containing individual ID (space) phenotype value");
	InitKey("WRITE_KNOWLEDGE_BIN_REPORT", "YES",  "Generates a text report describing the bins as contained within the knowledge layout");
	InitSingletary("BIN_COLLAPSE_THRESHOLD", "200",        "Set threshold for bin collapsing.");
	InitSingletary("WRITE_KNOWLEDGE_BINS", "YES",			   "Write bin output using knowledge trees");
	InitSingletary("INTERGENIC_BIN_LENGTH", "50",          "The number of kilobases intergenic bins can hold.");
	InitSingletary("BIN_TRAVERSE_THRESHOLD", "50",         "The minimum number of SNPs a group or region must contain before stepping down into a finer resolution.");
	InitSingletary("MINIMUM_BIN_SIZE",       "5",          "The minimum number of SNPs a bin must have to be analyzed.");
	InitSingletary("EXPAND_INTO_SUBGENES",   "YES",        "Should the traversal step into exons/introns?");
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
		tasks.insert(std::make_pair(item->taskType, item));
	else
		delete item;
}

void Configuration::ExecuteConfiguration(BinApplication* app) {
	//Write all of our settings to the relevant variables in memory
	//Knowledge::Region::DuplicateDD_Weight				= GetDouble("IMPLICATION_IDX_DUPLICATE_WEIGHT");
	//Knowledge::RegionManager::modelGenerationType	= Knowledge::ModelGenerationMode::ConvertType(GetString("DISEASE_DEPENDENT_LEVEL").c_str());
	//Knowledge::BinaryArchive								= GetBoolean("BINARY_MODEL_ARCHIVE");
	//Task::detailedReport				= GetBoolean("DETAILED_REPORTS");
	//Knowledge::GeneGeneModelArchive::minImplicationIndex		= GetInteger("MINIMUM_IMPLICATION_INDEX");
	//Knowledge::GeneGeneModelArchive::maxModelCount	= GetInteger("MAX_SNP_MODEL_COUNT");

	BinManager::mafCutoff									= GetDouble("MAF_CUTOFF");
	DataImporter::CompressedVCF							= GetBoolean("COMPRESSED_VCF");
	
	//Knowledge::KbGroup::CollapseAssociationReport = GetBoolean("COLLAPSE_ASSOCIATION_REPORT");
	app->SetGeneExtension(GetInteger("GENE_BOUNDARY_EXTENSION"));
	app->SetReportPrefix(GetString("REPORT_PREFIX").c_str());
	app->UseHtmlReports(GetBoolean("HTML_REPORTS"));
	BinManager::maxSnpCount = GetInteger("BIN_COLLAPSE_THRESHOLD");
	std::string sep = GetLine("OUTPUT_DELIMITER");
	if (sep.length() > 0) {
	    boost::replace_all(sep, std::string("\""), "");
	    boost::replace_all(sep, std::string("'"), "");
	}
	if (sep.length() == 0)
		sep = " ";
	SetValue("OUTPUT_DELIMITER", std::string("'")+sep+"'");
	Task::GenerateFiles::OutputDelimeter = sep;

	vector<string> pheno_files;
	GetLines("PHENOTYPE_FILENAME", pheno_files);

	Task::GenerateFiles::WriteBinData					= GetBoolean("WRITE_BIN_DATA");
	Task::GenerateFiles::WriteGenotypeData				= GetBoolean("WRITE_GENOTYPE_DATA");

	//TODO: Check to see what these options actually did
	//Task::BinCollapse::VisualizeGroupTrees				= GetBoolean("WRITE_COLLAPSABLE_BIN_REPORT");
	//Task::BinCollapse::WriteKnowledgeBins				= GetBoolean("WRITE_KNOWLEDGE_BINS");

	//Build out the task list
	if (Task::GenerateFiles::WriteBinData || Task::GenerateFiles::WriteGenotypeData) {
		Task::GenerateFiles *t = new Task::GenerateFiles();
		tasks.insert(TaskPair(t->taskType, t));
	}
	//AddTask("GENE_REPORT", new Biofilter::Task::GeneReport());
	//AddTask("MARKER_INFO_REPORT", new Biofilter::Task::MarkerInfo());
	//AddTask("SNP_GENE_MAP", new Biofilter::Task::SnpReport());
	//AddTask("SNP_REPORT", new Biofilter::Task::SnpGeneMap());
	//AddTask("EXPORT_GENE_MODELS", new Biofilter::Task::GeneGeneModelReport());
	//AddTask("EXPORT_SNP_MODELS", new Biofilter::Task::SnpSnpModelArchive());
	

	BinManager::IntergenicBinWidth = GetInteger("INTERGENIC_BIN_LENGTH")*1000;
	BinManager::BinTraverseThreshold = GetInteger("BIN_TRAVERSE_THRESHOLD");
	BinManager::MinBinSize = GetInteger("MINIMUM_BIN_SIZE");
	BinManager::ExpandByExons = GetInteger("EXPAND_INTO_SUBGENES");
	
	/*if (Task::BinCollapse::VisualizeGroupTrees || Task::BinCollapse::WriteKnowledgeBins) {
		Task::BinCollapse *item = new Task::BinCollapse();
		tasks.insert(TaskPair(item->taskType, item));
	}*/
	
	TaskList::iterator itr = tasks.begin();
	TaskList::iterator end = tasks.end();

	while (itr != end)
		itr++->second->Init(app);

	
}

}
