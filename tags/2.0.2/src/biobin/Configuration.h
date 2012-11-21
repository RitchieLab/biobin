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

namespace BioBin{

class Configuration{

public:
	class Bool{
	public:
		Bool() : _data(false) {}
		Bool(const bool d) : _data(d){}
		Bool(const std::string& d) : _data(fromString(d)){}

		operator const char*() const{ return _data ? "Y" : "N"; }
		operator bool() const{ return _data;}

	private:
		static bool fromString(const std::string& s){
			return (s[0] == 'Y' || s[0] == 'y');
		}

		bool _data;
	};

	static boost::program_options::options_description& addCmdLine(boost::program_options::options_description& opts);
	static boost::program_options::options_description& addConfigFile(boost::program_options::options_description& opts);
	static boost::program_options::options_description& addVisible(boost::program_options::options_description& opts);

	static void printConfig(std::ostream& os);
	static void printHelp(std::ostream& os);

	static void parseOptions(const boost::program_options::variables_map& vm);

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
	static void validate(boost::any& v, const std::vector<std::string>& values, std::vector<T>*, int);

	// Options allowed in both command line and config file
	static boost::program_options::options_description _generic;
	// Options allowed only in the config file
	static boost::program_options::options_description _config;
	// Options allowed only in the command line
	static boost::program_options::options_description _cmd;
	// Options allowed in config file and command line, but hidden from help
	static boost::program_options::options_description _hidden;

	static bool _generic_init;
	static bool _config_init;
	static bool _cmd_init;
	static bool _hidden_init;

};


template <class T>
void Configuration::validate(boost::any& v, const std::vector<std::string>& values, std::vector<T>*, int){
	if(v.empty()){
		v = boost::any(std::vector<T>());
	}

	std::vector<T>* val_vec_ptr = boost::any_cast<std::vector<T> >(&v);

	std::vector<std::string>::const_iterator itr = values.begin();
	std::vector<std::string>::const_iterator end = values.end();
	while(itr != end){
		val_vec_ptr->push_back(boost::lexical_cast<T>(*itr));
		++itr;
	}
}

}

namespace std{
ostream& operator<<(ostream& o, const BioBin::Configuration::Bool& d);
istream& operator>>(istream& in, BioBin::Configuration::Bool& d_out);
}


#endif /* CONFIGURATION_H_ */
