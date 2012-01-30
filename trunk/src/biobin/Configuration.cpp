/*
 * Configuration.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: jrw32
 */

#include "Configuration.h"
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>

#include "main.h"
#include "binmanager.h"
#include "taskfilegeneration.h"
#include "application.h"
#include "PopulationManager.h"

using std::string;
using std::vector;
using po::value;

using boost::shared_ptr;

namespace BioBin{

po::options_description Configuration::_generic("BioBin Options");
po::options_description Configuration::_cmd;
po::options_description Configuration::_hidden;
po::options_description Configuration::_config;

bool Configuration::_generic_init = false;
bool Configuration::_config_init = false;
bool Configuration::_cmd_init = false;
bool Configuration::_hidden_init = false;

void Configuration::initGeneric(){
	_generic.add_options()
		("settings-db,D", value<string>(&Main::c_knowledge_file)->default_value("knowledge.bio"),
				"The location of the database")
		("vcf-file,V",value<string>(&Main::c_vcf_file), "The file containing VCF information")
		("compressed-vcf,C", "Flag indicating VCF file is compressed")
		("maf-cutoff,F",value<float>(&BinManager::mafCutoff)->default_value(.05),
				"The maximum minor allele frequency to consider eligible for bin inclusion")
		("include-sources",value<vector<int> >()->composing()->multitoken(),
				"A list of source IDs to include")
		("include-source-file",value<vector<string> >()->composing()->multitoken(),
				"A list of filenames containing source IDs to include")
		("add-group", value<vector<string> >()->composing(),
				"A list of filenames containing a group collection definition")
		("output-delimiter,d",value<string>(&Task::GenerateFiles::OutputDelimiter)->default_value(","),
				"The delimiter to use when outputting text files")
		("phenotype-filename,p",value<vector<string> >(&PopulationManager::c_phenotype_files)->composing(),
				"Filename containing phenotype information")
		("bin-minimum-size,m", value<uint>(&BinManager::MinBinSize)->default_value(5),
				"The minimum size of any bin")
		("bin-expand-size,e", value<uint>(&BinManager::BinTraverseThreshold)->default_value(50),
				"The size above which bins are expanded into child bins, if possible")
		("bin-expand-exons,x","Flag indicating to expand bins into exons/intron regions")
		("bin-expand-functions,f","Flag indicating to expand bins by functionality of the variants")
		("intergenic-bin-length,i",value<uint>(&BinManager::IntergenicBinWidth)->default_value(50),
				"Number of kilobases intergenic bins can hold")
		("no-loci,l","Flag indicating desire to not write locus report")
		("no-bins,b","Flag indicating desire to not write bin report")
		("no-genotypes,g","Flag indicating desire to not write genotype report")
		("no-locus-freq,q","Flag indicating desire to not write Case v. Control Minor Allele Freq. report")
		("genomic-build,G",value<string>(&Main::c_genome_build)->default_value("37"),
				"Genomic build of input data")
		("phenotype-control-value", value<float>(&PopulationManager::c_phenotype_control),
				"Phenotype control value")
		("min-control-frac", value<float>(&PopulationManager::c_min_control_frac)->default_value(0.125),
				"Minimum fraction of population needed for control cases");

	_generic_init = true;
}

void Configuration::initHidden(){
	_hidden_init = true;
}

void Configuration::initCmd(){
	_cmd_init = true;
}

void Configuration::initConfig(){
	_config_init = true;
}

void Configuration::initAll(){
	if (!_cmd_init){
		initCmd();
	}
	if(!_generic_init){
		initGeneric();
	}
	if(!_hidden_init){
		initHidden();
	}
	if(!_config_init){
		initConfig();
	}
}

po::options_description& Configuration::addCmdLine(po::options_description& opts){
	initAll();
	return opts.add(_generic).add(_cmd).add(_hidden);
}

po::options_description& Configuration::addConfigFile(po::options_description& opts){
	initAll();
	return opts.add(_generic).add(_config).add(_hidden);
}

void Configuration::printConfig(std::ostream& os){
	po::options_description conf_options("Biobin Configuration Options");
	addConfigFile(conf_options);

	os << "############################\n";
	os << "# Biobin configuration options\n";
	os << "############################\n\n";

	vector<shared_ptr<po::option_description> >::const_iterator itr =
			conf_options.options().begin();
	vector<shared_ptr<po::option_description> >::const_iterator end =
				conf_options.options().end();

	boost::any dummy;
	// String to hold the text value of the default argument
	string txt_val;
	int eq_pos;

	while(itr != end){
		os << "# " << (*itr)->description() << std::endl;
		txt_val = "";
		if ((*itr)->semantic()->apply_default(dummy)){
			//If here, we have a default argument, so let's find out what it is!
			string name = (*itr)->semantic()->name();
			eq_pos = name.find_last_of('=');
			if(name.find_last_of('=') != string::npos){
				eq_pos = name.find_last_of('=');
				txt_val = name.substr(eq_pos + 1, name.find_last_of(')') - eq_pos - 1);
			}
		}else{
			os << "#";
		}

		os << (*itr)->long_name();
		if (txt_val.size()){
			os << " = " << txt_val;
		}

		os << std::endl << std::endl;

		++itr;
	}
}

void Configuration::parseOptions(const po::variables_map& vm){
	// NOTE: Help, version and sample will be parsed elsewhere!
	// In fact, they'll probably be moved out of here to where they're parsed

	if(vm.count("bin-expand-exons")){
		BinManager::ExpandByExons = true;
		std::cerr<<"WARNING: expansion into exons is not currently supported.\n";
	}

	if(vm.count("bin-expand-functions")){
		BinManager::ExpandByFunction = true;
		std::cerr<<"WARNING: expansion by functionality is not currently supported.\n";
	}

	if(vm.count("add-groups")){
		Main::c_custom_groups = vm["add-groups"].as<vector<string> >();
	}

	if(vm.count("vcf-file")){
		string fn(boost::filesystem::path(vm["vcf-file"].as<string>()).filename().string());
		Application::reportPrefix = fn.substr(0,fn.find_first_of('.'));
	}

	if(vm.count("include-sources")){
		std::cerr<<"WARNING: include-sources functionality has not been implemented yet.\n";
	}

	if(vm.count("include-source-file")){
		std::cerr<<"WARNING: include-source-file functionality has not been implemented yet.\n";
	}

	if(vm.count("no-loci")){
		BioBin::Task::GenerateFiles::WriteLociData = false;
	}

	if(vm.count("no-bins")){
		BioBin::Task::GenerateFiles::WriteBinData = false;
	}

	if(vm.count("no-genotypes")){
		BioBin::Task::GenerateFiles::WriteGenotypeData = false;
	}

	if(vm.count("no-locus-freq")){
		BioBin::Task::GenerateFiles::WriteAFData = false;
	}

}

}





