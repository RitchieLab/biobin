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
		/*
		std::ofstream file(filename.c_str());
		Utility::StringArray binNames = app->BinNames();
		file << "ID" << OutputDelimeter << "Status" << OutputDelimeter
				<< Utility::Join(binNames, OutputDelimeter.c_str()) << "\n";
		std::vector<uint> binHits;
		app->GetMaxBinHits(binHits); ///< This is number of SNPs
		std::vector<uint>::iterator bhitr = binHits.begin();
		std::vector<uint>::iterator bhend = binHits.end();

		while (bhitr != bhend) {
			*bhitr = *bhitr *= 2;
			bhitr++;
		}
		file << "Totals" << OutputDelimeter << uniqueStatus.size()
				<< OutputDelimeter
				<< Utility::Join(binHits, OutputDelimeter.c_str()) << "\n";
		while (itr != end) {
			itr->WriteBins(file, OutputDelimeter.c_str());
			file << "\n";
			itr++;
		}
		*/
	}
	if (WriteGenotypeData) {
		std::string filename = app->AddReport("data", "genotypes",
				"Genotype Data");
		app->writeGenotypeData(filename);
		/*
		std::ofstream file(filename.c_str());
		itr = individuals.begin();
		while (itr != end) {
			itr->WriteGenotypes(file, OutputDelimeter.c_str());
			file << "\n";
			itr++;
		}
		*/
	}
}

}
}
