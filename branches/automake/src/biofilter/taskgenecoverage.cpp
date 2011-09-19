/* 
 * File:   genecoverage.cpp
 * Author: torstees
 * 
 * Created on April 14, 2011, 11:58 AM
 */

#include "taskgenecoverage.h"

namespace Biofilter {
	namespace Task {

void GeneCoverage::GenerateTxtReport(const char *filename) {
	if (datasets.size() > 0) {
		if (Task::detailedReport)
			return GenerateDetailedTxtReport(filename);
		std::ofstream file(filename);
		uint regionCount = regions->Size();
		file<<"Gene,Total";
		std::map<std::string, Knowledge::SnpDataset>::iterator itr = datasets.begin();
		std::map<std::string, Knowledge::SnpDataset>::iterator end = datasets.end();

		while (itr != end) {
			file<<","<<itr++->first;
		}
		file<<"\n";

		for (uint i=0; i<regionCount; i++) {
			std::vector<int> counts;
			Knowledge::Region &r = (*regions)[i];
			std::set<Utility::Locus> totalSnps;

			snps->RangeSnpLookup(r.chrom, r.effStart, r.effEnd, totalSnps);

			counts.push_back(totalSnps.size());

			itr = datasets.begin();

			while (itr != end) {
				std::set<Utility::Locus> localSnps;
				itr->second.RangeSnpLookup(r.chrom, r.effStart, r.effEnd, localSnps);
				Utility::SnpArray common;

				std::set_intersection(totalSnps.begin(),
						totalSnps.end(),
						localSnps.begin(),
						localSnps.end(),
						std::inserter(common, common.begin()));

				counts.push_back(common.size());
				itr++;
			}
			file<<r.name<<","<<Utility::Join<std::vector<int> >(counts, ",")<<"\n";
		}
	} else {
		std::cerr<<"Coverage Reports require at least 1 coverage platform\n";
	}
}
void GeneCoverage::GenerateDetailedTxtReport(const char *filename) {
	std::ofstream file(filename);
	uint regionCount = regions->Size();
	file<<"Gene,Total";
	std::map<std::string, Knowledge::SnpDataset>::iterator itr = datasets.begin();
	std::map<std::string, Knowledge::SnpDataset>::iterator end = datasets.end();

	while (itr != end) {
		file<<","<<itr->first;
		file<<","<<itr++->first<<" SNPs";
	}
	file<<"\n";

	for (uint i=0; i<regionCount; i++) {
		Knowledge::Region &r = (*regions)[i];
		std::set<Utility::Locus> totalSnps;

		snps->RangeSnpLookup(r.chrom, r.effStart, r.effEnd, totalSnps);

		file<<r.name<<","<<totalSnps.size();
		std::set<Utility::Locus>::iterator sItr = totalSnps.begin();
		std::set<Utility::Locus>::iterator sEnd = totalSnps.end();
		uint count = 0;
		file<<",";
		while (sItr != sEnd)  {
			if (count++)
				file<<":";
			file<<sItr++->RSID();
		}


		itr = datasets.begin();

		while (itr != end) {
			std::set<Utility::Locus> localSnps;
			itr->second.RangeSnpLookup(r.chrom, r.effStart, r.effEnd, localSnps);
			std::set<Utility::Locus> common;

			std::set_intersection(totalSnps.begin(),
					totalSnps.end(),
					localSnps.begin(),
					localSnps.end(),
					std::inserter(common, common.begin()));

			std::set<Utility::Locus>::iterator cItr = common.begin();
			std::set<Utility::Locus>::iterator cEnd = common.end();
			uint ccount = 0;
			file<<","<<common.size()<<",";
			while (cItr != cEnd)  {
				if (ccount++)
					file<<":";
				file<<cItr++->RSID();
			}
			itr++;
		}
		file<<"\n";

	}
}


void GeneCoverage::AddSources(Utility::StringArray& rsFiles, Utility::StringArray& mapFiles, LiftOver::Converter* cnv) {
	if (rsFiles.size() + mapFiles.size() > 0) {
		std::set<std::string> snpsLost;			// With mapfiles, this becomes inadequate-probably should switch to a string set

		Utility::StringArray::iterator itr = rsFiles.begin();
		Utility::StringArray::iterator end = rsFiles.end();

		std::set<std::string> rsIDs;
		std::cerr<<"Loading Source for Gene Coverage:\n";
		while (itr != end) {
			std::string filename = *itr++;
			std::string fileContents = Utility::LoadContents(filename.c_str());
			Utility::StringArray lines = Utility::RsIDs(fileContents.c_str(), "\n");
			
			rsIDs.insert(lines.begin(), lines.end());
			std::cerr<<"* Platform        "<<filename<<": "<<rsIDs.size()<<" SNPs Loaded\n";
		}

	
		std::set<std::string> snpsFound;
		if (rsFiles.size() > 0) {
			snps->LoadData(rsIDs, snpsFound);
			std::set_difference(rsIDs.begin(), rsIDs.end(), snpsFound.begin(), snpsFound.end(), std::inserter(snpsLost, snpsLost.begin()));
		}
		
		//OK, now we create the RS datasets. We should do this before processing
		//map files, since we don't want to pollute the RS Maps with something 
		//pulled from the map file (which probably won't matter...but just in case)
		itr = rsFiles.begin();
		end = rsFiles.end();

		while (itr != end) {
			std::string filename				= *itr++;
			std::string fileContents		= Utility::LoadContents(filename.c_str());
			Utility::StringArray lines = Utility::RsIDs(fileContents.c_str(), "\n");
			std::set<std::string> rsIDs;
			rsIDs.insert(lines.begin(), lines.end());
			snps->GetFragment(rsIDs, datasets[Utility::ExtractBaseFilename(filename.c_str())]);
		}

		//Map files
		itr = mapFiles.begin();
		end = mapFiles.end();
		while (itr != end) {
			std::string filename = *itr++;
			Knowledge::SnpDataset &dataset = datasets[filename];
			std::multimap<Utility::Locus, Utility::Locus> converted;
			
			cnv->ConvertDataset(filename.c_str(), converted);
			std::multimap<Utility::Locus, Utility::Locus>::iterator citr = converted.begin();
			std::multimap<Utility::Locus, Utility::Locus>::iterator cend = converted.end();
			
			while (citr != cend) {
				if (citr->second.pos == 0)
					snpsLost.insert(citr->first.RSID());
				else {
					dataset.AddSNP(citr->second.chrom, citr->second.pos, citr->second.RSID().c_str());
					snps->AddSNP(citr->second.chrom, citr->second.pos, citr->second.RSID().c_str());
				}
				citr++;
			}			
		}


		std::cerr<<"Total SNPs found: "<<snps->Size()<<"\n";
		std::cerr<<"SNPs Lost       : "<<snpsLost.size()<<"\n";

	} else {
		throw Utility::Exception::General("Gene Coverage files require at least one RS or MAP based platform. Please see the manual for details.\n");
	}
}

}
}

