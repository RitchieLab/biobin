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
bool GenerateFiles::WriteAFData = true;
bool GenerateFiles::WriteBinFreqData = true;
std::string GenerateFiles::OutputDelimiter = ",";

GenerateFiles::GenerateFiles(BinApplication* app) :
		Task(3, app), _app(app) {
}

void GenerateFiles::ExecuteTask() {
	std::cout << "Executing Dataset File Generation\n";

	if (WriteBinData) {
		std::string filename = _app->AddReport("bins", "csv", "Bin Counts");
		_app->writeBinData(filename,OutputDelimiter);
	}
	if (WriteGenotypeData) {
		std::string filename = _app->AddReport("genotypes", "csv",
				"Genotype Data");
		_app->writeGenotypeData(filename,OutputDelimiter);
	}
	if (WriteLociData){
		std::string filename = _app->AddReport("locus", "csv",
				"Locus Data");
		_app->writeLoci(filename,OutputDelimiter);
	}
	if (WriteAFData){
		std::string filename = _app->AddReport("AllFreq", "csv",
				"Case vs. Control Allele Freq.");
		_app->writeAFData(filename,OutputDelimiter);
	}
	if(WriteBinFreqData){
		std::string filename = _app->AddReport("BinFreq", "csv",
						"Case vs. Control Bin Freq.");
		_app->writeBinFreqData(filename,OutputDelimiter);
	}
}

}
}
