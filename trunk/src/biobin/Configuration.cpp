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
using boost::any;

using namespace boost::program_options;

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
		("compressed-vcf,C", value<Bool>()->default_value(false), "Flag indicating VCF file is compressed")
		("maf-cutoff,F",value<float>(&BinManager::mafCutoff)->default_value(.05),
				"The maximum minor allele frequency to consider eligible for bin inclusion")
		("keep-common-loci,k",value<Bool>()->default_value(true),
				"Flag indicating to keep data pertaining to common variants (turn off to save memory)")
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
		("bin-expand-exons,x",value<Bool>()->default_value(false),
				"Flag indicating to expand bins into exons/intron regions")
		("bin-expand-functions,f",value<Bool>()->default_value(false),
				"Flag indicating to expand bins by functionality of the variants")
		("bin-pathways",value<Bool>()->default_value(true),
				"Flag indicating not to include pathways in the analysis")
		("bin-genes",value<Bool>()->default_value(true),
				"Flag indicating not to include genes in the analysis")
		("bin-intergenic",value<Bool>()->default_value(true),
				"Flag indicating not to include intergenic bins in the analysis")
		("intergenic-bin-length,i",value<uint>(&BinManager::IntergenicBinWidth)->default_value(50),
				"Number of kilobases intergenic bins can hold")
		("report-prefix",value<string>(), "A prefix to give to all of the reports")
		("report-loci",value<Bool>()->default_value(true),
				"Flag indicating desire to write locus report")
		("report-bins",value<Bool>()->default_value(true),
				"Flag indicating desire to write bin report")
		("report-genotypes",value<Bool>()->default_value(false),
				"Flag indicating desire to write genotype report")
		("report-locus-freq",value<Bool>()->default_value(false),
				"Flag indicating desire to write Case v. Control Minor Allele Freq. report")
		("report-bin-freq",value<Bool>()->default_value(false),
				"Flag indicating desire to write Bin Case v. Control Frequency report")
		("genomic-build,G",value<string>(&Main::c_genome_build)->default_value("37"),
				"Genomic build of input data")
		("phenotype-control-value", value<float>(&PopulationManager::c_phenotype_control),
				"Phenotype control value")
		("min-control-frac", value<float>(&PopulationManager::c_min_control_frac)->default_value(0.125),
				"Minimum fraction of population needed for control cases")
		("rare-case-control", value<Bool>()->default_value(false),
				"Flag indicating determining rarity of variants by both case and control populations")
		("disease-model",value<PopulationManager::DiseaseModel>(&PopulationManager::c_model)->default_value(PopulationManager::ADDITIVE),
				"Disease model (additive, dominant, or recessive)");

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


	if(vm.count("vcf-file")){
		string fn(boost::filesystem::path(vm["vcf-file"].as<string>()).filename().string());
		if(vm.count("report-prefix")){
			Application::reportPrefix = vm["report-prefix"].as<string>();
		}else{
			Application::reportPrefix = fn.substr(0,fn.find_first_of('.'));
		}
	}
	DataImporter::CompressedVCF = vm["compressed-vcf"].as<Bool>();
	DataImporter::KeepCommonLoci = vm["keep-common-loci"].as<Bool>();
	DataImporter::RareCaseControl = vm["rare-case-control"].as<Bool>();


	if(vm.count("add-groups")){
		Main::c_custom_groups = vm["add-groups"].as<vector<string> >();
	}

	if(vm.count("include-sources")){
		std::cerr<<"WARNING: include-sources functionality has not been implemented yet.\n";
	}

	if(vm.count("include-source-file")){
		std::cerr<<"WARNING: include-source-file functionality has not been implemented yet.\n";
	}

	//==========================================
	// Parsing report strategies
	//==========================================
	BioBin::Task::GenerateFiles::WriteLociData = vm["report-loci"].as<Bool>();
	BioBin::Task::GenerateFiles::WriteBinData = vm["report-bins"].as<Bool>();
	BioBin::Task::GenerateFiles::WriteGenotypeData = vm["report-genotypes"].as<Bool>();
	BioBin::Task::GenerateFiles::WriteAFData = vm["report-locus-freq"].as<Bool>();
	BioBin::Task::GenerateFiles::WriteBinFreqData = vm["report-bin-freq"].as<Bool>();

	//===========================================
	// Parsing binning strategies
	//===========================================
	BinManager::ExpandByExons = vm["bin-expand-exons"].as<Bool>();
	if(BinManager::ExpandByExons){
		std::cerr<<"WARNING: expansion into exons is not currently supported.\n";
	}
	BinManager::ExpandByFunction = vm["bin-expand-functions"].as<Bool>();
	if(BinManager::ExpandByFunction){
		std::cerr<<"WARNING: expansion by functionality is not currently supported.\n";
	}
	BioBin::BinManager::UsePathways = vm["bin-pathways"].as<Bool>();
	BioBin::BinManager::ExpandByGenes = vm["bin-genes"].as<Bool>();
	BioBin::BinManager::IncludeIntergenic = vm["bin-intergenic"].as<Bool>();
	if(!BioBin::BinManager::UsePathways && !BioBin::BinManager::ExpandByGenes){
		if(!BioBin::BinManager::IncludeIntergenic){
			std::cerr << "ERROR: You must bin by either pathways, genes, or intergenic.\n";
			throw validation_error(validation_error::invalid_option_value);
		}else{
			std::cerr<<"WARNING: You elected not to bin by pathways or genes.  " <<
							"You will only get intergenic bins.\n";
		}


	}
}

}


std::istream& operator>>(std::istream& in, BioBin::PopulationManager::DiseaseModel& model_out)
{
    std::string token;
    in >> token;
    if(token.size() > 0){
    	char s = token[0];
    	if(s == 'a' || s == 'A'){
    		model_out = BioBin::PopulationManager::ADDITIVE;
    	}else if(s == 'd' || s == 'D'){
    		model_out = BioBin::PopulationManager::DOMINANT;
    	}else if(s == 'r' || s == 'R'){
    		model_out = BioBin::PopulationManager::RECESSIVE;
    	}else{
    		throw validation_error(validation_error::invalid_option_value);
    	}
    }else{
    	throw validation_error(validation_error::invalid_option_value);
    }
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

namespace std{
std::istream& operator>>(std::istream& in, BioBin::Configuration::Bool& d_out)
{
    std::string token;
    in >> token;
    if(token.size() > 0){
    	char s = token[0];
    	d_out = (s == 'y' || s == 'Y');
    }else{
    	throw validation_error(validation_error::invalid_option_value);
    }
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

ostream& operator<<(ostream& o, const BioBin::Configuration::Bool& d){
	o << (const char *) d;
	return o;
}
}

