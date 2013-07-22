/* 
 * File:   ldspline.cpp
 * Author: torstees
 * 
 * Created on August 30, 2010, 12:16 PM
 */

#include "ldspline.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "timestamp.h"

using namespace Spline;


int main(int argc, char **argv) {
	if (argc < 3) {
		std::cerr<<"ldspline "<<APPMAJOR<<"."<<APPMINOR<<"."<<APPBUGFIX<<" ("<<BUILD_NUMBER<<") "<<BUILD_TYPE<<"  "<<BUILD_DATE<<"\n\n";
		std::cerr<<"Usage: ldspline [load/list/report/summarize/report-bounds/range/expand-bounds] filename [dp/rs] [min-ld_value] [chrom] [position|filename]..\n";
		std::cerr<<"\tload         - initiates loading LD information in haploview format.\n";
		std::cerr<<"\tlist         - lists the spline results for a given set of positions. \n";
		std::cerr<<"\treport-bounds- lists spline boundaries for a give set of positions.\n";
		std::cerr<<"\tsummarize    - lists basic details about loci contained within the data store. Summarize requires a chromosome number of ALL\n";
		std::cerr<<"\texpand-bounds- expands boundaries based on splines contained within the range\n";
		std::cerr<<"\t\tThe first parameter after list should be the file containing the LD data being queried.\n";
		std::cerr<<"\tdp/rs        - (list only, required) Indicates which type of statistic is being queried for (required only for list)\n";
		std::cerr<<"\tmin-ld_value - (list only, required) Indicates the threshold for the spline.\n";
		std::cerr<<"\tchrom        - (list only, required) Indicates which chromosome the positions are from.\n";
		std::cerr<<"\n* When running lists, users supply positions instead of filenames\n";
		return 1;
	}

	std::string cmd = argv[1];
	std::string dataFilename = LDUtility::StripExtension(argv[2]);
	std::string rawFilename = argv[2];

	if (cmd == "load") {
		for (int i=2; i<argc; i++) {
			LdSpline ldspline;
			std::string filename = ldspline.ConvertHaploLD(rawFilename.c_str());
			//ldspline.RunReport(std::cerr);
		}
	}



	else if (cmd == "list") {
		if (argc > 6) {
			LdSpline ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string type = argv[3];
			float ldvalue = atof(argv[4]);
			std::string chrom = argv[5];
			std::cout<<"chrom\tprimary rs\tprimary pos\tsecondary rs\tsecondary pos\t"<<type<<"/tdirection\n";
			for (int i=6; i<argc; i++) {
				ldspline.RunReport(chrom.c_str(), atoi(argv[i]), ldvalue, type.c_str(), std::cout);
			}
		}
		else {
			std::cerr<<"-- Insufficient parameters for list\n";
		}
	}
	else if (cmd == "summarize") {
		LdSpline ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
		std::string chrom = "ALL";
		if (argc > 3)
			chrom = argv[3];
		ldspline.Summarize(std::cerr, chrom.c_str());
	}
	else if (cmd == "expand-bounds") {
		if (argc > 6) {
			LdSpline ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string type = argv[3];
			float ldvalue = atof(argv[4]);
			std::string chrom = argv[5];
			std::cout <<"chrom\tinit_lower\tinit_upper\tlower_bound\tupper_bound\t"<<type<<"\n";
			for (int i=6; i<argc; i+=2) {
				if (atoi(argv[i]) < atoi(argv[i+1])) {
					if (type == "RS" || type == "rs") {
						std::pair<int, int> bounds = ldspline.GetBoundariesRS(chrom.c_str(), atoi(argv[i]), atoi(argv[i+1]), ldvalue);
						std::cerr<<chrom<<"\t"<<argv[i]<<"\t"<<argv[i+1]<<"\t"<<bounds.first<<"\t"<<bounds.second<<"\t"<<ldvalue<<"\n";
					}
					else if (type == "DP" || type == "dp") {
						std::pair<int, int> bounds = ldspline.GetBoundariesDP(chrom.c_str(), atoi(argv[i]), atoi(argv[i+1]), ldvalue);
						std::cerr<<chrom<<"\t"<<argv[i]<<"\t"<<argv[i+1]<<"\t"<<bounds.first<<"\t"<<bounds.second<<"\t"<<ldvalue<<"\n";
					}
				} else {
					std::cerr<<"\nLower bound must actually be lower than the upper bound: "<<argv[i]<<" > "<<argv[i+1]<<"\n";
					return 1;
				}
			}
		}
	}
	else if (cmd == "report-bounds") {
		if (argc > 6) {
			LdSpline ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string type = argv[3];
			float ldvalue = atof(argv[4]);
			std::string chrom = argv[5];
			std::cout<<"chrom\trs\tpos\tlower_pos\tlower_rs\tupper_bound\tupper_rs\t"<<type<<"\n";
			for (int i=6; i<argc; i++) {
				ldspline.ReportBounds(chrom.c_str(), atoi(argv[i]), ldvalue, type.c_str(), std::cout);
			}
		}
		else {
			std::cerr<<"-- Insufficient parameters for list\n";
		}

	}
	else if (cmd == "report") {
		LdSpline ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
		ldspline.RunReport(std::cerr);
	}
	else if (cmd == "export-lomap") {
		LdSpline ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str(), true);
		std::string loFilename = LDUtility::ExtractBaseFilename(argv[2]) + ".bim";
		ldspline.ExportForLiftOver(loFilename.c_str());
		std::cerr<<"Lift over output: "<<loFilename<<"\n";
	}
	else if (cmd == "import-lomap") {
		LdSpline ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str(), true);

		std::string loFilename = LDUtility::ExtractBaseFilename(argv[2]) + ".new";
		std::string loUnmapped = LDUtility::ExtractBaseFilename(argv[2]) + ".unmapped";
		std::cerr<<"Loading position data from "<<loFilename<<"\t"<<loUnmapped<<"\n";
		ldspline.ImportLiftOver(loFilename.c_str(), loUnmapped.c_str());
		ldspline.SaveToCopyBinary(std::string(dataFilename + "-b" + std::string(argv[3]) + ".ldspline").c_str());
	}
	else if (cmd == "asdf") {
		for (int i=2; i<argc; i++) {
			LdSpline ldspline2;
			LdSpline ldspline;
			std::string filename = ldspline.ConvertHaploLD(argv[i]);
			ldspline2.LoadFromBinary(filename.c_str());
			ldspline.RunReport(std::cerr);
			ldspline2.RunReport(std::cerr);
		}
	}
	else {
		std::cerr<<"Unknown command: "<<argv[1]<<" Options include: load|list|summarize\n";
	}




	return 0;
}
