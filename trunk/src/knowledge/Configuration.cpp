/*
 * Configuration.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: jrw32
 */

#include "Configuration.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#include "RegionCollection.h"
#include "GroupCollection.h"

using std::string;
using std::vector;
using po::value;
using boost::shared_ptr;
using std::ifstream;
using boost::algorithm::split;
using boost::algorithm::is_any_of;

namespace Knowledge{

po::options_description Configuration::_generic;
po::options_description Configuration::_cmd;
po::options_description Configuration::_hidden;
po::options_description Configuration::_config;

bool Configuration::_generic_init = false;
bool Configuration::_config_init = false;
bool Configuration::_cmd_init = false;
bool Configuration::_hidden_init = false;

void Configuration::initGeneric(){
	_generic.add_options()
			("include-groups",value<vector<int> >()->composing()->multitoken(),
					"A list of group IDs to include")
			("include-group-names", value<vector<string> >()->composing()->multitoken(),
					"A list of group names to include")
			("include-group-file", value<vector<string> >()->composing(),
					"A file containg a group definition")
			("population,P", value<string>(&RegionCollection::pop_str)->default_value("NO-LD"),
					"The population to base the gene boundaries on")
			("gene-boundary-extension,B", value<int>(&RegionCollection::gene_expansion)->default_value(0),
					"The amount to expand the genes by (when using NO-LD)");

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
	po::options_description conf_options("LOKI Configuration Options");
	addConfigFile(conf_options);

	os << "############################\n";
	os << "# LOKI configuration options\n";
	os << "############################\n\n";

	vector<shared_ptr<po::option_description> >::const_iterator itr =
			conf_options.options().begin();
	vector<shared_ptr<po::option_description> >::const_iterator end =
				conf_options.options().end();

	boost::any dummy;
	// String to hold the text value of the default argument
	string txt_val;
	while(itr != end){
		os << "# " << (*itr)->description() << std::endl;
		txt_val = "";
		if ((*itr)->semantic()->apply_default(dummy)){
			//If here, we have a default argument, so let's find out what it is!
			string name = (*itr)->semantic()->name();
			if(name.find_last_of('=') != string::npos){
				int eq_pos = name.find_last_of('=');
				txt_val = name.substr(eq_pos - 1, name.find_last_of(')') - eq_pos - 1);
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
	if (vm.count("include-group-names")){
		GroupCollection::c_group_names = vm["include-group-names"].as<vector<string> >();
	}

	if (vm.count("include-groups")){
		vector<int> ids = vm["include-groups"].as<vector<int> >();
		GroupCollection::c_id_list.insert(ids.begin(), ids.end());
	}

	if (vm.count("include-group-file")){
		vector<string> files = vm["include-group-file"].as<vector<string> >();
		// Iterate over each file, adding it to the group names
		vector<string>::const_iterator f_itr = files.begin();
		vector<string>::const_iterator f_end = files.end();
		while(f_itr != f_end){
			ifstream data_file((*f_itr).c_str());
			if (!data_file.is_open()){
				std::cerr<<"WARNING: cannot find " << (*f_itr) <<", ignoring.";
			} else {
				string line;
				vector<string> result;
				while(data_file.good()){
					getline(data_file, line);
					split(result, line, is_any_of(" \n\t"));
					GroupCollection::c_group_names.insert(GroupCollection::c_group_names.end(), result.begin(), result.end());
				}
				data_file.close();
			}
			++f_itr;
		}
	}
}
}



