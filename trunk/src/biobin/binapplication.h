/* 
 * File:   binapplication.h
 * Author: torstees
 *
 * Created on June 22, 2011, 10:35 AM
 * 
 * One of the big TODO list things would be to integrate the two 
 * forms of SNPs: Biofilter and BioBin. Right now, we have two
 * very different approaches to SNPs. For the biofilter, we only
 * need a way to recognize names and associate them with a basepair
 * and chromosome. However, for biobin, we need to maintain alleles
 * and provide the ability to perform genotype conversion and 
 * some other stuff. So, the Locus object is much more complex. 
 * 
 * Most likely it's just a matter of moving the biobin Locus class
 * to someplace common and changing the biofilter code to use it...
 * 
 */

#ifndef BINAPPLICATION_H
#define	BINAPPLICATION_H

#include "application.h"

#include "Bin.h"

#include "binmanager.h"
#include "PopulationManager.h"

#include <utility>
#include <string>
#include <vector>
#include <set>
#include <map>

using std::string;
using std::vector;
using std::map;
using std::set;

#include "knowledge/Region.h"
#include "knowledge/Locus.h"
#include "knowledge/Group.h"
#include "knowledge/liftover/ConverterSQLite.h"

namespace BioBin {
	
class BinApplication : public Application {
public:
	BinApplication(const string& vcf_file);
	virtual ~BinApplication() {}

	template <class SNP_cont>
	void InitVcfDataset(string& genomicBuild,
			SNP_cont& lostSnps);

	/**
	 * Initialize the bins.  After this call, the binmanager will have the final
	 * bins set up and ready to output.
     */
	void InitBins();

	void writeBinData(const string& filename, const string& sep=",") const;
	void writeGenotypeData(const string& filename, const string& sep=",") const;
	void writeLoci(const string& filename, const string& sep=",") const;
	void writeAFData(const string& filename, const string& sep=",") const;
	void writeBinFreqData(const string& filename, const string& sep=",") const;



	
private:
	PopulationManager _pop_mgr;

	BinManager binData;

};


template <class SNP_cont>
void BinApplication::InitVcfDataset(std::string& genomicBuild, SNP_cont& lostSnps) {

	// First things first, let's load our individuals
	//const vector<bool>& controls = _pop_mgr.loadIndividuals(_data);

	//locusRemap.clear();
	//locusRemap.reserve(locusArray.size());
	Knowledge::Liftover::ConverterSQLite cnv(genomicBuild, _db);
	int chainCount = cnv.Load();

	if (chainCount > 0) {
		_pop_mgr.loadLoci(dataset,&cnv);
	}else{
		_pop_mgr.loadLoci(dataset);
	}
	/*

		vector<Knowledge::Locus*> toDelete;
		vector <Knowledge::Locus*> locusArray;
		//map <Knowledge::Locus*, vector<short> > tmp_genotype_map;
		_data.getLoci(locusArray, controls);

		std::string conversionLog = this->AddReport("lift-over", "tsv", "SNPs that were lifted over to new build which differed dramatically or changed chromosome");
		std::ofstream cnvLog(conversionLog.c_str());
		cnvLog<<"Chrom(Orig),Pos(Orig),RSID(Old),Chrom(New),Pos(New),RSID(New)\n";

		map<Knowledge::Locus*, Knowledge::Locus*> converted;
		vector<Knowledge::Locus*> not_found;
		cnv.convertLoci(locusArray.begin(), locusArray.end(), converted, not_found);
		vector<Knowledge::Locus*>::iterator itr = locusArray.begin();
		vector<Knowledge::Locus*>::iterator end = locusArray.end();
		map<Knowledge::Locus*, Knowledge::Locus*>::const_iterator missing =
				converted.end();
		map<Knowledge::Locus*, Knowledge::Locus*>::const_iterator new_loc_itr;

		std::stringstream missingSNPs;
		while (itr != end) {

			Knowledge::Locus &orig = **itr;
			new_loc_itr = converted.find(&orig);

			if (new_loc_itr != missing) {

				if(!((*new_loc_itr).second)){
					missingSNPs << *((*new_loc_itr).first) << "\n";
				} else {
					if ((*new_loc_itr).first->getChrom() != (*new_loc_itr).second->getChrom()) {
						cnvLog << *((*new_loc_itr).first) << ","
								<< *((*new_loc_itr).second) << "\n";
					}

					Locus* newLoc = (*new_loc_itr).second;
					dataset.push_back(newLoc);
					_data.remapLocus(&orig, newLoc);
				}
			}else{
				missingSNPs << orig << "\n";
			}

			++itr;
		}
		if (missingSNPs.str().length() > 0) {
			std::string filename = AddReport("missing-snps", "txt", "SNPs that were dropped during build conversion");
			std::ofstream file(filename.c_str());
			file<<missingSNPs.str();
		}

		for (vector<Knowledge::Locus*>::iterator del_it = locusArray.begin(); del_it != locusArray.end(); del_it++){
			delete *del_it;
		}

	} else {
		_data.getLoci(dataset, controls);
	}

	_pop_mgr.loadGenotypes(dataset, _data);
	*/

}

}
#endif	/* BINAPPLICATION_H */

