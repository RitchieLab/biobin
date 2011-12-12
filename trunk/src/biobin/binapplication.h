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
//#include "individual.h"

#include "dataimporter.h"

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
	BinApplication();
	virtual ~BinApplication();

	template <class SNP_cont>
	void InitVcfDataset(string& genomicBuild,
			SNP_cont& lostSnps,
			DataImporter& vcfimporter);

	/**
	 * Initialize the bins.  After this call, the binmanager will have the final
	 * bins set up and ready to output.
     */
	void InitBins();

	//void ApplyPhenotypes();

	template<class T_cont>
	void loadPhenotypes(const T_cont& pheno_files){_pop_mgr.loadPhenotypes(pheno_files);}

	void writeBinData(const string& filename, const string& sep=",") const;
	void writeGenotypeData(const string& filename, const string& sep=",") const;
	void writeLoci(const string& filename, const string& sep=",") const;


	/**
	 * Return the individuals that have been loaded 
    * @return
    */
	//const vector<Individual>& Individuals();

	/**
	 * Returns the number of SNPs that might contribute to a given bin
    * @param hits
    */
	//void GetMaxBinHits(vector<uint>& hits);

	//void GetBinContributors(vector<set<uint> >& contributors);

	/**
	 * Returns an array matching each of the bin names
    * @return
    */
	//const Utility::StringArray& BinNames();

	/**
	 * Returns the region object at a given index
    * @param idx
    * @return
    */
	//const Knowledge::Region& GetRegion(uint idx);


	//Utility::Locus &Locus(uint idx);

	/**
	 * Passes a list of genotypes (for all people at a given SNP) and returns the bin IDs 
    * @param snpIndex which SNP we are referring to
    * @param genotypes the original data from the vcf files
    * @param hits Vector containing the individual indexes where a variation occured
    * @param genotypes This is where we'll write genotype data
    * @return set of bin indexes for which this SNP applies
    */
	//std::set<uint> ParseSNP(uint snpIndex, vector<char>& genotypes, vector<Individual>& data);
	

	/**
	 * returns lookup region index -> snp index
	 */
	//void GenerateBinContentLookup(multimap<uint, uint>& binContents);
	//map<uint, uint> GetBinLookup();
private:
	///< Used to build and parse data into bins and genotypes
	// Also holds the actual bins!
	BinManager binData;

	PopulationManager _pop_mgr;


						///< Help to extract genotype data from vcf files
	//vector<Individual> individuals;	///< This represents the actual data from the vcf files
	//set<uint> binnable;					///< Just a map into the regions that created the entries
	//map<uint, uint> binIndex;


	//std::vector<Utility::Locus> loci;		///< The loci associated with the dataset
	//Knowledge::SnpDataset loci;			/// This is now dataset
};


template <class SNP_cont>
void BinApplication::InitVcfDataset(std::string& genomicBuild, SNP_cont& lostSnps, DataImporter& vcfimporter) {

	std::cerr<<"Loading VCF Data\n";

	//locusRemap.clear();
	//locusRemap.reserve(locusArray.size());
	Knowledge::Liftover::ConverterSQLite cnv(genomicBuild, _db);
	int chainCount = cnv.Load();

	if (chainCount > 0) {
		vector <Knowledge::Locus*> locusArray;
		//map <Knowledge::Locus*, vector<short> > tmp_genotype_map;
		vcfimporter.getLoci(locusArray);

		std::string conversionLog = this->AddReport("lift-over", "tsv", "SNPs that were lifted over to new build which differed dramatically or changed chromosome");
		std::ofstream cnvLog(conversionLog.c_str());
		cnvLog<<"Chrom(Orig),Pos(Orig),RSID(Old),Chrom(New),Pos(New),RSID(New)\n";

		multimap<Knowledge::Locus*, Knowledge::Locus*> converted;
		cnv.ConvertDataset(locusArray.begin(), locusArray.end(), converted);
		vector<Knowledge::Locus*>::iterator itr = locusArray.begin();
		vector<Knowledge::Locus*>::iterator end = locusArray.end();
		multimap<Knowledge::Locus*, Knowledge::Locus*>::const_iterator missing =
				converted.end();


		uint i=0;
		//uint validLocus = 0;
		std::stringstream missingSNPs;
		while (itr != end) {
			Knowledge::Locus &orig = **itr;

			if (converted.find(&orig) != missing) {
				multimap<Knowledge::Locus*, Knowledge::Locus*>::const_iterator map_itr =
						converted.lower_bound(&orig);
				multimap<Knowledge::Locus*, Knowledge::Locus*>::const_iterator last  =
						converted.upper_bound(&orig);

				if (converted.count(&orig) != 1){
					std::cerr<<"It was observed that there are multiple hits returned by convert dataset: "<<orig.getID()<<" has "<<converted.count(&orig)<<" counterparts.\n";
				}

				while (map_itr != last) {
					if (!((*map_itr).second) || (*map_itr).second->getPos() == 0){
						missingSNPs<<*((*map_itr).first)<<"\n";
					}else {
						if ((*map_itr).second->getChrom() != -1) {
							if ((*map_itr).first->getChrom() != (*map_itr).second->getChrom() ||
									((float)abs((float)(*map_itr).first->getPos() - (float)(*map_itr).second->getPos())/(float)(*map_itr).first->getPos())> 0.01){
								cnvLog<<*((*map_itr).first)<<","
										<<*((*map_itr).second)<<"\n";
							}

							dataset.push_back((*map_itr).second);
							//dataset.AddSNP(first->second.chrom, itr->second.pos, first->second.RSID().c_str());
							//locusRemap[itr->second.chrom].push_back(validLocus++);
							//locusArray[i] = itr->second;
						}
						else{
							cnvLog<<*((*map_itr).first)<<","
							<<" **** \n";
						}
					}
					++map_itr;
				}
			}
			delete *itr;
			++itr;
			++i;
		}
		if (missingSNPs.str().length() > 0) {
			std::string filename = AddReport("missing-snps", "txt", "SNPs that were dropped during build conversion");
			std::ofstream file(filename.c_str());
			file<<missingSNPs.str();
		}
	} else {
		vcfimporter.getLoci(dataset);
	}

	// OK, now that we've loaded the dataset, we can load up the individuals
	_pop_mgr.loadIndividuals(vcfimporter);
	_pop_mgr.loadGenotypes(dataset, vcfimporter);

}

}
#endif	/* BINAPPLICATION_H */

