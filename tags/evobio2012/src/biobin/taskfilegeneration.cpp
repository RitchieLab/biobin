/* 
 * File:   taskfilegeneration.cpp
 * Author: torstees
 * 
 * Created on June 30, 2011, 9:50 AM
 */

#include "taskfilegeneration.h"

#include "binapplication.h"

namespace BioBin {
namespace Task {

bool GenerateFiles::WriteBinData = true;
bool GenerateFiles::WriteGenotypeData = true;
bool GenerateFiles::WriteLociData = true;
std::string GenerateFiles::OutputDelimiter = ",";

GenerateFiles::GenerateFiles(BinApplication* app) :
		Task(3, app), _app(app) {
}

void GenerateFiles::ExecuteTask() {
	std::cerr << "Executing Dataset File Generation\n";

	if (WriteBinData) {
		std::string filename = _app->AddReport("data", "bins", "Bin Counts");
		_app->writeBinData(filename);
	}
	if (WriteGenotypeData) {
		std::string filename = _app->AddReport("data", "genotypes",
				"Genotype Data");
		_app->writeGenotypeData(filename);
	}
	if (WriteLociData){
		std::string filename = _app->AddReport("locus", "csv",
				"Locus Data");
		_app->writeLoci(filename);
	}
}

}
}
