/* 
 * File:   binmanager.h
 * Author: torstees
 *
 * Created on July 21, 2011, 3:54 PM
 */

#ifndef BINMANAGER_H
#define	BINMANAGER_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>

//#include "individual.h"
#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Group.h"
#include "knowledge/Locus.h"

#include "Bin.h"

using std::map;
using std::vector;
using std::string;
using std::pair;
using std::set;

namespace BioBin {

class PopulationManager;

class BinManager {
public:
	typedef set<Bin*>::const_iterator const_iterator;

	BinManager(const PopulationManager& pop_mgr);

	virtual ~BinManager();
	//BinManager(const BinManager& orig);
	
	void InitBins(const map<uint, Knowledge::GroupCollection*> &groups,
			const Knowledge::RegionCollection& regions,
			const vector<Knowledge::Locus*>& loci);

	int numRareVariants() const { return _rare_variants.size();}
	int numVariants() const {return _total_variants;}
	int numBins() const {return _bin_list.size();}
	//int numTotalVariants() const {return _total_variants;}

	const_iterator begin() const {return _bin_list.begin();}
	const_iterator end() const {return _bin_list.end();}

	void printBins(std::ostream& os, Knowledge::Locus* locus, const string& sep=":") const;

	static uint IntergenicBinWidth;				///< The width of the intergenic bins within a chromosome
	static uint BinTraverseThreshold;			///< The number of SNPs to determine whether we continue traversing
	static uint MinBinSize;							///< How small do we tolerate bins to be before we ignore the bin altogether
	static bool ExpandByGenes;						///< Do we want to drop down to genes, if the group is large enough?
	static bool ExpandByExons;						///< Do we want to drop to to introns and exons, if the group is large enough
	static bool ExpandByFunction;					///< Indicate that we do want to use function 
	static float mafCutoff; 	///< Max maf to produce result in a bin

	// Carried this over from taskbincollapse, but I'm not sure what it does
	static uint maxSnpCount;
private:

	BinManager(const BinManager&);
	BinManager& operator=(const BinManager&);

	// Collapses all of the bins according to the preferences we set.
	void collapseBins();
	// Erases (and increments) a single bin
	void eraseBin(set<Bin*>::iterator&);

	Bin* addRegionBin(Knowledge::Region* reg);

	// The authoritative list of all of the bins.  Everything else holds pointers
	// to bins in this set.
	set<Bin*> _bin_list;
	// Mappint of Region IDs to bins
	map<int, Bin*> _region_bins;
	// Mapping of group IDs to bins
	map<int, Bin*> _group_bins;
	// List of intergenic bins
	map<pair<short, int>, Bin*> _intergenic_bins;
	// List of bins by locus
	map<Knowledge::Locus*, set<Bin*> > _locus_bins;

	// List of all rare variants
	set<Knowledge::Locus*> _rare_variants;


	int _total_variants;

	const PopulationManager& _pop_mgr;


};


} //namespace BioBin

#endif	/* BINMANAGER_H */

