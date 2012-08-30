/*
 * Configuration.h
 *
 *  Created on: Dec 14, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_CONFIGURATION_H
#define KNOWLEDGE_CONFIGURATION_H

#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace po = boost::program_options;

using std::stringstream;
using std::vector;
using std::string;
using std::map;
using boost::tokenizer;
using boost::escaped_list_separator;
using po::validation_error;

namespace Knowledge{

class Configuration{

public:

	template <class T>
	class Container{
	public:
		explicit Container<T>(const string& str){stringstream ss(str);ss >> (*this);}
		Container<T>(){}
		operator vector<T>() const {return _data;}
		void push_back(T val){ _data.push_back(val);}
	private:
		vector<T> _data;
	};


	static po::options_description& addCmdLine(po::options_description& opts);
	static po::options_description& addConfigFile(po::options_description& opts);
	static po::options_description& addVisible(po::options_description& opts);

	static void printConfig(std::ostream& os);
	//static void printHelp(std::ostream& os);
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

}

namespace std{

template <class T>
ostream& operator<<(ostream& o, const Knowledge::Configuration::Container<T>& d){
	string sep = ",";
	const vector<T>& data = (vector<T>) d ;
	//vector<T>::const_iterator itr = data.begin();
	for(unsigned int i=0; i<data.size(); i++){
		if(i){
			o << sep;
		}
		o << data[i];
	}
	o << "\n";

	return o;
}

template <class T>
istream& operator>>(istream& in, Knowledge::Configuration::Container<T>& d_out){
	string in_s;
	string sep = ",";
	in >> in_s;
	// Now, split up the string
	tokenizer<escaped_list_separator<char> > tok(in_s);
	for(tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
		try{
			 d_out.push_back(boost::lexical_cast<T>(*beg));
		}catch(boost::bad_lexical_cast& e){
	    	throw validation_error(validation_error::invalid_option_value);
		}
	}

	return in;
}

}


#endif /* CONFIGURATION_H_ */
