/* 
 * File:   taskfilegeneration.cpp
 * Author: torstees
 * 
 * Created on June 30, 2011, 9:50 AM
 */

#include "taskfilegeneration.h"
#include "utility/strings.h"
namespace BioBin {
namespace Task {

bool GenerateFiles::WriteBinData = true;
bool GenerateFiles::WriteGenotypeData = true;
std::string GenerateFiles::OutputDelimeter = " ";

GenerateFiles::GenerateFiles() :
		Task(3), app(NULL) {
}

GenerateFiles::~GenerateFiles() {
}

void GenerateFiles::Init(Application* app) {
	this->app = (BinApplication*) app;
}

void GenerateFiles::ExecuteTask() {
	std::cerr << "Executing Dataset File Generation\n";

	if (WriteBinData) {
		std::string filename = app->AddReport("data", "bins", "Bin Counts");
		app->writeBinData(filename);
	}
	if (WriteGenotypeData) {
		std::string filename = app->AddReport("data", "genotypes",
				"Genotype Data");
		app->writeGenotypeData(filename);
	}
	std::string filename = app->AddReport("locus", "csv",
			"Locus Data");
	app->writeLoci(filename);
}

}
}
