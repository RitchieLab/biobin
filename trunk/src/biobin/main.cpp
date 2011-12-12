#include "main.h"
#include <iostream>
#include "utility/exception.h"

/**
 * The VCF Tools require this, but I don't want to use their main...
 */
namespace VCF {
	std::ofstream LOG;
}

namespace BioBin {


void Main::LoadConfiguration(const char *cfgFilename) {
	cfg.SetValue("REPORT_PREFIX", Utility::ExtractBaseFilename(cfgFilename));
	cfg.Parse(cfgFilename);
}


void Main::InitRegionData() {
	vector<string> missingAliases;
	vector<string> aliasList;

	app.LoadRegionData(cfg.GetLine("POPULATION"), missingAliases, aliasList);
}

// TODO: deal w/ schizophrenic vector/string uint/int decisions
void Main::InitGroupData() {
	//User defined groups
	Utility::StringArray udGroups;
	cfg.GetLines("ADD_GROUP", udGroups);

	//Any specialized searches are defined here
	vector<string> lines;
	vector<uint> groupIDs;
	cfg.GetIntegers("INCLUDE_GROUPS", groupIDs);
//	std::cerr<<"ASDFASDF: "<<lines.size()<<"\n";
	set<uint> ids;
	ids.insert(groupIDs.begin(), groupIDs.end());

	cfg.LoadFileContents("INCLUDE_GROUP_FILE", ids);

	//Now, let's do the same for names
	
	vector<string> groups;
	cfg.GetLines("INCLUDE_GROUP_NAMES", groups);
	cfg.LoadFileContents("INCLUDE_GROUP_NAME_FILE", groups);

	vector<int> group_id_input;
	group_id_input.reserve(ids.size());
	set<uint>::const_iterator itr = ids.begin();
	set<uint>::const_iterator end = ids.end();
	while(itr != end){
		group_id_input.push_back(*itr);
		++itr;
	}
	app.LoadGroupDataByName(udGroups, groups, group_id_input);

}

void Main::RunCommands() {
	VCF::LOG.open("vcf-responses.log");

	app.Init(cfg.GetLine("SETTINGS_DB"), !silentRun);
	std::string genomicBuild = cfg.GetString("GENOMIC_BUILD");
	if (genomicBuild != "") {
		app.LoadBuildConverter(genomicBuild);
	}
	/* This is biofilter specific, not for biobin
	switch (action) {
		case BiofilterAction::ListGroups:
			{
				std::vector<std::string> keywords;
				std::string s = cfg.GetLine("GROUP_SEARCH_CRITERIA");
				if (s != "All")
					keywords = Utility::Split(s.c_str(), ",");
				
				app.ListGroupIDs(std::cout, keywords);
			} return;
		case BiofilterAction::ListPopulationIDs:
			app.ListPopulationIDs(std::cout);
			return;
		case BiofilterAction::ListGenes: {
			std::string s = cfg.GetLine("GENE_COVERAGE");
			Utility::StringArray aliasList;
			if (s != "ALL")
				aliasList = Utility::Split(s.c_str(), ",");
			Utility::StringArray aliasTypeList;
			s = cfg.GetLine("ALIAS_TYPES");
			if (s != "ALL")
				aliasTypeList = Utility::Split(s.c_str(), ",");

			app.ListGenes(std::cout, aliasList, aliasTypeList);
			return;
		}
		case BiofilterAction::ListMetaGroups:{
//				InitGroupData();
//				app.ListMetaGroups(cout);
//
			return;
		}
	}*/
	
	//Tasks that run before SNPs load (not sure what those would be)
	cfg.RunTasks(0);

	vector<string> genes;
	std::string geneFilename = cfg.GetLine("GENE_COVERAGE");
	if (geneFilename != "")
		cfg.LoadFileContents("GENE_COVERAGE", genes);
	/**
	 * Do the SNP oriented stuff here
    */

	//With some new decisions about binning, we have to 
	//preload the groups-then we will handle the dataset merging/squeeze
	//This sort of makes tasks type 2 and 3 redundant....and is backward
	//from our previous approach

	// TODO: check for existence of the file here!
	string fn = cfg.GetLine("VCF_FILE");
	if(fn.size() == 0){
		std::cerr<<"No SNP dataset available. Unable to continue.\n";
		exit(1);
	}

	DataImporter vcfimporter(fn);
	LoadSNPs(vcfimporter);

	cfg.RunTasks(1);

	InitRegionData();

	cfg.RunTasks(2);

	InitGroupData();
	app.InitBins();

	cfg.RunTasks(3);

	//We need to make sure that there is one or more tasks at level four before we generate the models
	// We should NOT be here in biobin!!  There's no model generation going on!
	//if (cfg.CountTasks(4)) {
	//	app.ProduceModels(std::cout);
	//	cfg.RunTasks(4);
	//}
}



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

void Main::PrintBanner()  {
	if (!silentRun) {
		//std::cerr<<"biobin "<<APPMAJOR<<"."<<APPMINOR<<"."
		//	 <<APPBUGFIX<<" ("<<BUILD_NUMBER<<") "<<BUILD_TYPE<<"  "<<BUILD_DATE<<"\n";
#ifdef USE_MPI
		cout<<"* This application is compiled to run on parallel computing systems using MPI\n";
##else
		cout<<"* (serial)\n";
#endif
		std::cerr<<"\nMarylyn Ritchie, Carrie Buchanan and Eric Torstenson\n\n";
	}
}

void Main::PrintHelp() {
	silentRun = false;
	PrintBanner();
#ifdef USE_MPI
	cerr<<"usage: biofilter <configuration file> [ [command] ...] [ [parameter] ...]\n";
#else
	std::cerr<<"usage: biobin <configuration file> \n";
#endif
	std::cerr<<"\biobin - TODO write stuff about biobin\n";
	std::cerr<<"Optional Commands Include:\n";
	std::cerr<<"\t-S [--sample-config]                       -- Print sample configuration to std-out\n";
	std::cerr<<"\t--list-genes                               -- Lists all genes that are covered by at least one SNP\n";
	std::cerr<<"\t--genes alias_list alias_type              -- Lists all genes present in the database that match one of the comma \n"
		 <<"\t                                              separated. Either or both can also be ALL, which will show them all. \n";
	std::cerr<<"\nOptional Parameters Include:\n";
	std::cerr<<"\t-V [--vcf-file] <filename>                 -- Specify a vcf file to be used.\n";
	std::cerr<<"\t-B [--build] <build version>               -- Define the build associated with map files (35, 36, 37)\n";
	std::cerr<<"\t-P [--list-populations]                    -- Lists all available Population based LD boundary options\n";
	std::cerr<<"\t-p [--set-population] pop                  -- Override the configurations population setting (NO-LD, CEUDP1.0, etc)\n";
	std::cerr<<"\t-t [--maf-threshold] thresh                -- MAF cutoff for calling rare-variants.\n";
	std::cerr<<"\t-k [--knowledge-threshold] thresh          -- Max number of SNPs for a knowledge group to be binned.\n";
	std::cerr<<"\t--fix-variations var-filename-path         -- Sets the path (and filename) to the appropriate variation file.\n";
	std::cerr<<"\t                                              This should only be done if the file needs to be moved to a new location.\n";
}


void Main::LoadSNPs(DataImporter& vcf) {
		std::string genomicBuild = cfg.GetLine("GENOMIC_BUILD");
		vector<string> lostSnps;
		app.InitVcfDataset(genomicBuild, lostSnps, vcf);
		std::cerr<<lostSnps.size()<<" SNPs were not able to be found in the variations database.\n";
}

void Main::InitGroups() {

}

int Main::SetConfigValue(int nextCmd, int argc, const char *var, const char *val, const char *err) {
	if (nextCmd < argc) 
		cfg.SetValue(var, val);
	else {
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
	if (strcmp(argv[curr], "--DB")==0)
		return SetConfigValue(nextCmd, argc, "SETTINGS_DB", argv[nextCmd], "--DB must be followed by a database filename");
	if (strcmp(argv[curr], "--marker-info")==0) {
		cfg.SetValue("MARKER_INFO_REPORT", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "-b")==0 || strcmp(argv[curr], "--binary")==0)
		return SetConfigValue(nextCmd, argc, "BINARY_MODEL_ARCHIVE", argv[nextCmd], "--binary must be followed by Yes/No");
	if (strcmp(argv[curr], "-P")==0 || strcmp(argv[curr], "--list-populations")==0) {
		action = BiofilterAction::ListPopulationIDs;
		return nextCmd;
	}
	if (strcmp(argv[curr], "-D")==0) {
		cfg.SetValue("DETAILED_REPORTS", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "--map-snps-to-gene")==0) {
		cfg.SetValue("SNP_GENE_MAP", "ON");
		if (nextCmd < argc) {
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
		}
		else {
			action = BiofilterAction::ParseError;
			std::cerr<<"--map-snps-to-gene must be followed by a filename containing a list of genes.\n";
			return -1;
		}
		return nextCmd;
	}
	if (strcmp(argv[curr], "--list-genes")==0)  {
		cfg.SetValue("GENE_REPORT", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "--snp-report")==0)  {
		cfg.SetValue("SNP_REPORT", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "--map-snps-to-gene")==0)  {
		cfg.SetValue("SNP_GENE_MAP", "ON");
		return nextCmd;
	}
	if (strcmp(argv[curr], "-G")==0 || strcmp(argv[curr], "--groups")==0) {
		if (argc > nextCmd) {
			silentRun = true;				// At this point, we don't care about the other stuff
			cfg.SetValue("LIST_GROUPS_FROM_DB", "ON");
			cfg.SetValue("GROUP_SEARCH_CRITERIA", argv[nextCmd++]);
			action = BiofilterAction::ListGroups;
			return -1;
		}
		else {
			action = BiofilterAction::ParseError;
			std::cerr<<"--groups must include search criterion or ALL (to list all groups).\n";
			return -1;
		}
	}
	if (strcmp(argv[curr], "--genes")==0) {
		if (argc > nextCmd + 1) {
			silentRun = true;				// At this point, we don't care about the other stuff
			cfg.SetValue("LIST_GENES_FROM_DB", "ON");
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
			cfg.SetValue("ALIAS_TYPES", argv[nextCmd++]);
			action = BiofilterAction::ListGenes;
			return -1;
		}
		else {
			action = BiofilterAction::ParseError;
			std::cerr<<"--genes must include genes (comma separated) and alias type (comma separated). Either can be replaced by ALL.\n";
			return -1;
		}

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

	if (strcmp(argv[curr], "-k")==0 || strcmp(argv[curr], "--knowledge-threshold")==0)
		return SetConfigValue(nextCmd, argc, "BIN_COLLAPSE_THRESHOLD", argv[nextCmd], "--knowledge-threshold must be followed by the max number of SNPs allowed in knowledge based bins.");
	if (strcmp(argv[curr], "-t")==0 || strcmp(argv[curr], "--maf-threshold")==0) 
		return SetConfigValue(nextCmd, argc, "MAF_CUTOFF", argv[nextCmd], "--MAF_CUTOFF must be followed by maf threshold value");
	if (strcmp(argv[curr], "--PREFIX")==0)
		return SetConfigValue(nextCmd, argc, "REPORT_PREFIX", argv[nextCmd], "--PREFIX must be followed by prefix to be prepended to the generated filenames");
	if (strcmp(argv[curr], "-p")==0 || strcmp(argv[curr], "--set-population")==0)
		return SetConfigValue(nextCmd, argc, "POPULATION", argv[nextCmd], "--set-population must be followed by name population you wish to use");
	if (strcmp(argv[curr], "--gene-boundary")==0)
		return SetConfigValue(nextCmd, argc, "GENE_BOUNDARY_EXTENSION", argv[nextCmd], "--gene-boundary must be followed by an integer describing the number of bases");
	if (strcmp(argv[curr], "-V")==0 || strcmp(argv[curr], "--vcf-file")==0)
		return SetConfigValue(nextCmd, argc, "VCF_FILE", argv[nextCmd], "--vcf-file must be followed by the name of a vcf file.");
	if (strcmp(argv[curr], "-z")==0 || strcmp(argv[curr], "--gzipped vcf")==0)
		return SetConfigValue(nextCmd, argc, "COMPRESSED_VCF", argv[nextCmd], "--gzipped vcf must be followed by YES/NO.");
	action = BiofilterAction::ParseError;
	std::cerr<<"Unrecognized parameter: "<<argv[curr]<<"\n";
	return -1;
}



}


int main(int argc, char *argv[]) {
	std::string cfgFilename;

	BioBin::Main *app = new BioBin::Main();					///<The application object


	if (!app->ParseCmdLine(argc, argv)) {
		delete app;
		exit(1);
	}
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
