#ifndef MAIN_H
#define MAIN_H

#include "configuration.h"

namespace BioBin {
class Main {
public:
	Main();
	~Main();
	/**
	 * @brief Pass the arguments to the application object
	 */
	bool ParseCmdLine(int argc, char **argv);
	int ParseCmd(int curr, int argc, char **argv);
	void PrintHelp();						///< Display usage details
	void PrintBanner();					///< Display details about the software
	void LoadConfiguration(const char *cfgFilename);
	void RunCommands();
	int SetConfigValue(int nextCmd, int argc, const char *var, const char *val, const char *err);
	void InitGroupData();
	void InitRegionData();
protected:
	void LoadSNPs(std::vector<uint>& remap, DataImporter& vcf);
	void InitGroups();

	struct BiofilterAction {
		enum Action {
			NoAction,						///< No particular action. This is the default
			ParseError,						///< Error loading configuration
			ListGroups,						///< List gorup IDs (optionally based on search)
			ListAliasTypes,				///< Lists the Alias types in the database
			ListGenes,						///< Just a DB dump of the genes in the database (based on one or more aliases/ALL and one or more type/ALL)
			ListGenesSimple,				///< Lists the genes covered by one or more SNPs in the SNP source
			ListMetaGroups,				///< Display only metagroups and their IDs (useful for group filtering)
			ListPopulationIDs,			///< Lists available population IDs (for ld expansion)
			PrintSampleConfig,			///< Write sample configuration file
			ListAssociations,				///< gene->group associations for a snp list
			FilterByGenes,					///< Filter a SNP list by presence inside a gene
			SetVariationFilename			///< Used to fix the variation filename
		};
	};


	
	Configuration cfg;					///< Configuration settings
	BinApplication app;					///< The application that does all of the work
	bool silentRun;						///< Used to silence the banner banter
	BiofilterAction::Action action;	///< Actual task to be performed


};
inline
Main::Main() : silentRun(false), action(BiofilterAction::NoAction) { }

inline
Main::~Main() { }
}


#endif
