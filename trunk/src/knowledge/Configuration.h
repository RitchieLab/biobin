/*
 * Configuration.h
 *
 *  Created on: Dec 14, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_CONFIGURATION_H
#define KNOWLEDGE_CONFIGURATION_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


namespace Knowledge{

class Configuration{

public:

	template <class T>
	class Container{
	public:
		explicit Container<T>(const std::string& str){std::stringstream ss(str);ss >> (*this);}
		Container<T>(){}
		operator std::vector<T>() const {return _data;}
		void push_back(const T& val){ _data.push_back(val);}
	private:
		std::vector<T> _data;
	};


	static boost::program_options::options_description& addCmdLine(boost::program_options::options_description& opts);
	static boost::program_options::options_description& addConfigFile(boost::program_options::options_description& opts);
	static boost::program_options::options_description& addVisible(boost::program_options::options_description& opts);

	static void printConfig(std::ostream& os);
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

/** Validates sequences. Allows multiple values per option occurrence
   and multiple occurrences. */
template<class T>
void validate(boost::any& v,
              const std::vector< std::string >& s,
              Knowledge::Configuration::Container<T>*,
              int)
{
    if (v.empty()) {
        v = boost::any(Knowledge::Configuration::Container<T>());
    }
    Knowledge::Configuration::Container<T>* tv = boost::any_cast< Knowledge::Configuration::Container<T> >(&v);
    assert(NULL != tv);
    for (unsigned i = 0; i < s.size(); ++i)
    {
    	std::stringstream ss(s[i]);
    	ss >> (*tv);
    }
}

template <class T>
std::ostream& operator<<(std::ostream& o, const Configuration::Container<T>& d){
	std::string sep = ",";
	const std::vector<T>& data = (std::vector<T>) d ;
	//vector<T>::const_iterator itr = data.begin();
	for(unsigned int i=0; i<data.size(); i++){
		if(i){
			o << sep;
		}
		o << data[i];
	}
	//o << "\n";

	return o;
}

template <class T>
std::istream& operator>>(std::istream& in, Configuration::Container<T>& d_out){
	std::string in_s;
	in >> in_s;
	// Now, split up the string
	boost::tokenizer<boost::escaped_list_separator<char> > tok(in_s);
	for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
		try{
			T obj(boost::lexical_cast<T>(*beg));
			d_out.push_back(obj);
		}catch(boost::bad_lexical_cast& e){
	    	throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
		}
	}

	return in;
}

}

#endif /* CONFIGURATION_H_ */
