/*
 * main.cpp
 *
 *  Created on: Apr 23, 2012
 *      Author: jrw32
 */

#include "ldsplineimporter.h"
#include <string>
#include <iostream>

using std::string;
using std::cerr;

void printHelp(){
	cerr << "usage: pop_loader [OPTIONS]\n\n";
	cerr << "pop_loader is a helper function for the LOKI database that will create\n";
	cerr << "population-specific boundaries in LOKI from an LDSpline file.\n\n";
	cerr << "Options Include:\n";
	cerr << "\t-h [--help]               Show this help message\n";
	cerr << "\t-d [--DB] <filename>      Give the location of the LOKI database \n"
		 << "\t                          (default: knowledge.bio)\n";
	cerr << "\t-c [--config] <filename>  Give the location of the Config file\n"
	     << "\t                          (default: ldspline.cfg)\n";
}

int main(int argc, char** argv){
	string db_fn = "knowledge.bio";
	string config_fn = "ldspline.cfg";

	// Parse the command line here
	int curr = 0;
	while (++curr < argc){
		string arg = argv[curr];
		if((arg == "--DB" || arg == "-d") && curr < argc){
			db_fn = argv[++curr];
		}
		if((arg == "--config" || arg == "-c" ) && curr < argc){
			config_fn = argv[++curr];
		}
		if(arg == "--help" || arg == "-h"){
			printHelp();
			exit(1);
		}
	}


	try{
		LdSplineImporter lds(config_fn, db_fn);
		lds.loadPops();
	}catch(std::runtime_error& e){
		std::cerr<<"Caught error " << e.what() << "\n";
		return 1;
	}
	return 0;
	//
}



