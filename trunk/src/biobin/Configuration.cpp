/*
 * Configuration.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: jrw32
 */

#include "Configuration.h"
#include <string>
#include <vector>
#include <sstream>
#include <set>

#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "main.h"
#include "binmanager.h"
#include "binapplication.h"
#include "PopulationManager.h"

#include "tests/TestFactory.h"
#include "tests/Test.h"

using std::string;
using std::vector;
using std::ostream;
using std::istream;
using std::stringstream;
using std::set;

using boost::shared_ptr;
using boost::any;

using BioBin::Test::TestFactory;
//using BioBin::Test::Test;

namespace po = boost::program_options;

using po::value;

//using namespace boost::program_options;

namespace BioBin{

po::options_description Configuration::_generic("BioBin Options");
po::options_description Configuration::_cmd("Command Line Options");
po::options_description Configuration::_hidden;
po::options_description Configuration::_config;

bool Configuration::_generic_init = false;
bool Configuration::_config_init = false;
bool Configuration::_cmd_init = false;
bool Configuration::_hidden_init = false;

void Configuration::initGeneric(){

	po::options_description binning_opts("Bin Generation Options");
	binning_opts.add_options()
		("maf-cutoff,F",value<float>(&BinManager::mafCutoff)->default_value(.05, "0.05"),
				"The maximum minor allele frequency to consider eligible for bin inclusion")
		("maf-threshold", value<float>(&BinManager::mafThreshold)->default_value(0, "0"),
				"The minimum minor allele frequency to consider eligible for bin inclusion")
		("keep-monomorphic", value<Bool>()->default_value(false),
				"Keep monomorphic markers (will not contribute to any bin)")
		("bin-minimum-size,m", value<uint>(&BinManager::MinBinSize)->default_value(5),
				"The minimum size of any bin")
		("bin-expand-size,e", value<uint>(&BinManager::BinTraverseThreshold)->default_value(50),
				"The size above which bins are expanded into child bins, if possible")
		("bin-expand-roles,x",value<Bool>()->default_value(false),
				"Flag indicating to expand bins into exons/intron regions")
		("filter-bin-role",value<Bool>()->default_value(false),
				"Flag indicating desire to filter by unknown role")
		("keep-unknown-role",value<Bool>()->default_value(false),
				"When filtering by bin role, if true, keep only unknown role bins, if false, drop unknown role bins")
		("bin-pathways",value<Bool>()->default_value(false),
				"Flag indicating not to include pathways in the analysis")
		("bin-regions",value<Bool>()->default_value(true),
				"Flag indicating not to include genes in the analysis")
		("bin-interregion",value<Bool>()->default_value(true),
				"Flag indicating not to include intergenic bins in the analysis")
		("interregion-bin-length,i",value<unsigned int>(&BinManager::IntergenicBinWidth)->default_value(50),
				"Number of kilobases intergenic bins can hold")
		("interregion-bin-step",value<unsigned int>(),
				"Sliding distance for intergenic bins, in kilobases (default = interregion-bin-length)")

				;

	po::options_description report_options("Report Generation Options");
	report_options.add_options()
		("report-prefix",value<string>(), "A prefix to give to all of the reports")
		("report-loci",value<Bool>()->default_value(true),
				"Flag indicating desire to write locus report")
		("report-bins",value<Bool>()->default_value(true),
				"Flag indicating desire to write bin report")
		("transpose-bins", value<Bool>()->default_value(false),
				"Transpose the Bin report (bins on rows)")
		("no-summary", value<Bool>()->default_value(false),
				"Suppress the summary information in a Bin report")
		("output-delimiter,d",value<string>(&Main::OutputDelimiter)->default_value(","),
				"The delimiter to use when outputting text files")

		;

	po::options_description pheno_opts("Phenotype Options");
	pheno_opts.add_options()
		("phenotype-filename,p",value<string>(&PopulationManager::c_phenotype_file),
				"Filename containing phenotype information")
		("min-control-frac", value<float>(&PopulationManager::c_min_control_frac)->default_value(0.125,"0.125"),
				"Minimum fraction of population needed for control cases")
		("rare-case-control", value<Bool>()->default_value(true),
				"Flag indicating determining rarity of variants by both case and control populations")
		("phenotype-control-value", value<float>(&PopulationManager::c_phenotype_control)->default_value(0, "0"),
				"Phenotype control value")
		;

	stringstream test_ss;
	Test::TestFactory::const_iterator tf_itr = Test::TestFactory::getFactory().begin();
	Test::TestFactory::const_iterator tf_end = Test::TestFactory::getFactory().end();
	test_ss << "Statistical tests to run on bins (" << tf_itr->first;
	while(++tf_itr != tf_end){
		test_ss << "," << tf_itr->first;
	}
	test_ss << ")";

	po::options_description test_opts("Testing Options");
	test_opts.add_options()
		("covariates", value<string>(&PopulationManager::c_covariate_file),
				"Filename containing covariate information")
		("weight-loci", value<Bool>()->default_value(false),
				"Add weights to the Locus")
		("weight-model",value<PopulationManager::WeightModel>(&PopulationManager::c_weight_type)->default_value(PopulationManager::MIN),
				"Method of determining weight for Loci (maximum, minimum, control, or overall)")
		("disease-model",value<PopulationManager::DiseaseModel>(&PopulationManager::c_model)->default_value(PopulationManager::ADDITIVE),
				"Disease model (additive, dominant, or recessive)")
		("test", value<vector<string> >()->composing(),
				test_ss.str().c_str())
		;


	_generic.add_options()
		("settings-db,D", value<string>(&Main::c_knowledge_file)->default_value("knowledge.bio"),
				"The location of the database")
		("vcf-file,V",value<string>(&Main::c_vcf_file), "The file containing VCF information")
		("threads,t", value<unsigned int>(&BinApplication::n_threads)->default_value(1),
				"Number of threads to use when PheWAS binning")
		("add-group", value<vector<string> >()->composing(),
				"A list of filenames containing a group collection definition")
		("genomic-build,G",value<string>(&Main::c_genome_build)->default_value("37"),
				"Genomic build of input data")
		;

	_generic.add(binning_opts);
	_generic.add(pheno_opts);
	_generic.add(report_options);
	_generic.add(test_opts);


	_generic_init = true;
}

void Configuration::initHidden(){
	_hidden_init = true;
}

void Configuration::initCmd(){
	_cmd.add_options()
		("print-populations", "Print populations available in LOKI")
		("print-sources", "Print the sources available in LOKI");
	_cmd_init = true;
}

void Configuration::initConfig(){
	_config_init = true;
}

void Configuration::initAll(){
	if(!_generic_init){
		initGeneric();
	}
	if (!_cmd_init){
		initCmd();
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
	return opts.add(_cmd).add(_generic).add(_hidden);
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

	vector<boost::shared_ptr<po::option_description> >::const_iterator itr =
			conf_options.options().begin();
	vector<boost::shared_ptr<po::option_description> >::const_iterator end =
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
			BinApplication::reportPrefix = vm["report-prefix"].as<string>();
		}else{
			BinApplication::reportPrefix = fn.substr(0,fn.find_first_of('.'));
		}
	}
	PopulationManager::RareCaseControl = vm["rare-case-control"].as<Bool>();
	PopulationManager::c_use_calc_weight = vm["weight-loci"].as<Bool>();
	PopulationManager::c_keep_monomorphic = vm["keep-monomorphic"].as<Bool>();

	if(vm.count("add-groups")){
		Main::c_custom_groups = vm["add-groups"].as<vector<string> >();
	}

	//==========================================
	// Parsing printing options
	//==========================================
	BinApplication::c_print_populations = vm.count("print-populations");
	BinApplication::c_print_sources = vm.count("print-sources");
	if(BinApplication::s_run_normal){
		BinApplication::s_run_normal = !(BinApplication::c_print_populations || BinApplication::c_print_sources);
	}

	//==========================================
	// Parsing report strategies
	//==========================================
	BioBin::Main::WriteLociData = vm["report-loci"].as<Bool>();
	BioBin::Main::WriteBinData = vm["report-bins"].as<Bool>();
	BinApplication::c_transpose_bins = vm["transpose-bins"].as<Bool>();
	PopulationManager::NoSummary = vm["no-summary"].as<Bool>();

	//===========================================
	// Parsing binning strategies
	//===========================================
	BinManager::ExpandByExons = vm["bin-expand-roles"].as<Bool>();
	BinManager::FilterByRole = vm["filter-bin-role"].as<Bool>();
	BinManager::KeepUnknown = vm["keep-unknown-role"].as<Bool>();
	BinManager::UsePathways = vm["bin-pathways"].as<Bool>();
	BinManager::ExpandByGenes = vm["bin-regions"].as<Bool>();
	BinManager::IncludeIntergenic = vm["bin-interregion"].as<Bool>();

	BinManager::IntergenicBinStep = vm.count("interregion-bin-step") ?
				vm["interregion-bin-step"].as<unsigned int>() :
				BinManager::IntergenicBinWidth;

	if(BinApplication::s_run_normal && !BioBin::BinManager::UsePathways && !BioBin::BinManager::ExpandByGenes){
		if(!BioBin::BinManager::IncludeIntergenic){
			std::cerr << "ERROR: You must bin by either pathways, regions, or interregion.\n";
			throw po::validation_error(po::validation_error::invalid_option_value);
		}else{
			std::cerr<<"WARNING: You elected not to bin by pathways or regions.  " <<
							"You will only get interregion bins.\n";
		}
	}

	//===========================================
	// Parsing test strategies
	//===========================================

	if(vm.count("test")){
		set<string> seen_tests;
		vector<string> test_strs = vm["test"].as<vector<string> >();
		vector<string> test_names;
		for(unsigned int i=0; i<test_strs.size(); i++){
			test_names.clear();
			boost::algorithm::split(test_names, test_strs[i],
					boost::is_any_of(","), boost::token_compress_on);
			for(unsigned int j=0; j<test_names.size(); j++){
				if(!seen_tests.count(test_names[j])){
					Test::Test* t = TestFactory::getFactory().Create(test_names[j]);
					if(!t){
						std::cerr << "WARNING: Test '" << test_names[j] << "' "
								<<"is not a recognized test, ignoring." << std::endl;
					} else {
						PopulationManager::c_tests.push_back(t);
						seen_tests.insert(test_names[j]);
					}
				}
			}
		}
	}
}

//}


//namespace std{

std::istream& operator>>(std::istream& in, BioBin::Configuration::Bool& d_out)
{
    std::string token;
    in >> token;
    if(token.size() > 0){
    	char s = token[0];
    	d_out = (s == 'y' || s == 'Y');
    }else{
    	throw po::validation_error(po::validation_error::invalid_option_value);
    }
//    else throw boost::program_options::validation_error("Invalid unit");
    return in;
}

ostream& operator<<(ostream& o, const BioBin::Configuration::Bool& d){
	o << (const char *) d;
	return o;
}
}

