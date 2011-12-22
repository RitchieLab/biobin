/*
 * Configuration.h
 *
 *  Created on: Dec 15, 2011
 *      Author: jrw32
 */

#ifndef BIOBIN_CONFIGURATION_H
#define BIOBIN_CONFIGURATION_H

#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include <iostream>
#include <vector>
#include <string>

using std::string;
using std::vector;

namespace po = boost::program_options;

namespace BioBin{

class Configuration{

public:
	static po::options_description& addCmdLine(po::options_description& opts);
	static po::options_description& addConfigFile(po::options_description& opts);
	static po::options_description& addVisible(po::options_description& opts);

	static void printConfig(std::ostream& os);
	static void printHelp(std::ostream& os);
	//static void printOptions(std::ostream& os, const po::variables_map& vm);

	static void parseOptions(const po::variables_map& vm);

private:
	// No construction or assignment of this class.  EVER!
	Configuration();
	Configuration(const Configuration&);
	Configuration& operator=(const Configuration&);

	// Initializers for static members
	static void initGeneric();
	static void initCmd();
	static void initHidden();
	static void initConfig();

	static void initAll();

	template <class T>
	static void validate(boost::any& v, const vector<string>& values, vector<T>*, int);

	// Options allowed in both command line and config file
	static po::options_description _generic;
	// Options allowed only in the config file
	static po::options_description _config;
	// Options allowed only in the command line
	static po::options_description _cmd;
	// Options allowed in config file and command line, but hidden from help
	static po::options_description _hidden;

	static bool _generic_init;
	static bool _config_init;
	static bool _cmd_init;
	static bool _hidden_init;

};


template <class T>
void Configuration::validate(boost::any& v, const vector<string>& values, vector<T>*, int){
	if(v.empty()){
		v = boost::any(vector<T>());
	}

	vector<T>* val_vec_ptr = boost::any_cast<vector<T> >(&v);

	vector<string>::const_iterator itr = values.begin();
	vector<string>::const_iterator end = values.end();
	while(itr != end){
		val_vec_ptr->push_back(boost::lexical_cast<T>(*itr));
		++itr;
	}
}

}


#endif /* CONFIGURATION_H_ */
