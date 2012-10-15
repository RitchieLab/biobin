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
	BinApplication(const string& db_fn, const string& vcf_file);
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

	static bool c_transpose_bins;
	
private:
	void printEscapedString(ostream& os, const string& toPrint, const string& toRepl, const string& replStr) const;
	string getEscapeString(const string& sep) const;

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

}

}
#endif	/* BINAPPLICATION_H */

