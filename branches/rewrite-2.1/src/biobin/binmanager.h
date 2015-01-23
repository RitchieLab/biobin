/* 
 * File:   binmanager.h
 * Author: torstees
 *
 * Created on July 21, 2011, 3:54 PM
 */

#ifndef BIOBIN_BINMANAGER_H
#define	BIOBIN_BINMANAGER_H

#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <deque>
#include <boost/unordered_map.hpp>

#define BOOST_IOSTREAMS_USE_DEPRECATED
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/operations.hpp>

#include <cstdio>

#include "knowledge/GroupCollection.h"
#include "knowledge/RegionCollection.h"
#include "knowledge/Group.h"
#include "knowledge/Locus.h"

#include "Bin.h"
#include "PopulationManager.h"

namespace Knowledge{
class Information;
}

namespace BioBin {

class PopulationManager;

class BinManager {
public:
	typedef std::set<Bin*>::const_iterator const_iterator;

	BinManager(const PopulationManager& pop_mgr,
			const Knowledge::RegionCollection& regions,
			const std::deque<Knowledge::Locus*>& loci,
			const Knowledge::Information& info,
			const PopulationManager::Phenotype& pheno);

	virtual ~BinManager();
	//BinManager(const BinManager& orig);
	
	void InitBins(const std::deque<Knowledge::Locus*>& loci);

	int numRareVariants() const { return _rare_variants;}
	int numBins() const {return _bin_list.size();}
	//int numTotalVariants() const {return _total_variants;}

	const_iterator begin() const {return _bin_list.begin();}
	const_iterator end() const {return _bin_list.end();}

	void printBinData(std::ostream& os, const std::string& sep, bool transpose= false) const;

	// create a temporary file (using tmpfile) and
	template <class L_cont>
	FILE* printLocusBins(const L_cont& loci, const std::string& sep=":") const;
	void printBins(std::ostream& os, Knowledge::Locus* locus, const std::string& sep=":") const;
	void printLocusBinCount(std::ostream& os, float pct=0.1) const;

	static unsigned int IntergenicBinWidth;				///< The width of the intergenic bins within a chromosome
	static unsigned int IntergenicBinStep;				///< The size of the step to take for intergenic sliding-window analysis
	static unsigned int BinTraverseThreshold;			///< The number of SNPs to determine whether we continue traversing
	static unsigned int MinBinSize;							///< How small do we tolerate bins to be before we ignore the bin altogether
	static bool UsePathways;						///< Do we want to use pathways in the analysis?
	static bool IncludeIntergenic;					///< Include Intergenic bins?
	static bool ExpandByGenes;						///< Do we want to drop down to genes, if the group is large enough?
	static bool ExpandByExons;						///< Do we want to drop to to introns and exons, if the group is large enough
	//! do we want to filter the unknown role bins?
	static bool FilterByRole;
	//! Keep or drop the unknown bins?
	static bool KeepUnknown;
	static float mafCutoff; 	///< Max maf to produce result in a bin
	static float mafThreshold; ///< Min maf to produce result in a bin

private:

	BinManager(const BinManager&);
	BinManager& operator=(const BinManager&);

	// Collapses all of the bins according to the preferences we set.
	void collapseBins();
	// Erases (and increments) a single bin
	void eraseBin(std::set<Bin*>::iterator&);

	Bin* addRegionBin(Knowledge::Region* reg);

	// The authoritative list of all of the bins.  Everything else holds pointers
	// to bins in this set.
	std::set<Bin*> _bin_list;
	// Mapping of Region IDs to bins
	boost::unordered_map<int, Bin*> _region_bins;
	// Mapping of group IDs to bins
	boost::unordered_map<int, Bin*> _group_bins;
	// List of intergenic bins
	boost::unordered_map<std::pair<short, int>, Bin*> _intergenic_bins;
	// List of bins by locus
	boost::unordered_map<Knowledge::Locus*, std::set<Bin*> > _locus_bins;

	// # of all variants, rare and common
	int _rare_variants;

	const PopulationManager& _pop_mgr;
	const Knowledge::RegionCollection& _regions;
	const Knowledge::Information& _info;

	const PopulationManager::Phenotype& _pheno;
};


template <class L_cont>
FILE* BinManager::printLocusBins(const L_cont& loci, const std::string& sep) const{
	FILE* tmpf = std::tmpfile();
	boost::iostreams::stream<boost::iostreams::file_descriptor> tmp_stream(boost::iostreams::file_descriptor(fileno(tmpf)),
					std::ios_base::binary | std::ios_base::in | std::ios_base::out);

	tmp_stream << _pop_mgr.getPhenotypeName(_pheno.getIndex());
	typename L_cont::const_iterator l_itr = loci.begin();
	while(l_itr != loci.end()){
		tmp_stream << "\n";
		printBins(tmp_stream, *l_itr, sep);
		++l_itr;
	}

	return tmpf;
}

} //namespace BioBin

#endif	/* BINMANAGER_H */

