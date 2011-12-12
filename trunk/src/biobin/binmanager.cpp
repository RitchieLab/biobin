/* 
 * File:   binmanager.cpp
 * Author: torstees
 * 
 * Created on July 21, 2011, 3:54 PM
 */

#include "binmanager.h"

#include "Bin.h"

#include "knowledge/RegionCollection.h"
#include "knowledge/Region.h"

using Knowledge::RegionCollection;
using Knowledge::Region;
using std::make_pair;

namespace BioBin {

uint BinManager::IntergenicBinWidth = 50000;
uint BinManager::BinTraverseThreshold = 50;
uint BinManager::MinBinSize = 1;
bool BinManager::ExpandByExons = true;
bool BinManager::ExpandByFunction = true;
float BinManager::mafCutoff = 0.05;
uint BinManager::maxSnpCount = 200;

BinManager::BinManager() : _total_variants(0){
}

BinManager::~BinManager() {
	set<Bin*>::const_iterator itr = _bin_list.begin();
	set<Bin*>::const_iterator end = _bin_list.end();
	while(itr != end){
		delete *itr;
		++itr;
	}
}

/**
 * Basically, we'll need to generate a list of "leaf" nodes for
 * the selected set of groups. We'll store these as pointers to 
 * the group in a set or a vector.
 * 
 * If the user wishes to drill down further (to the gene level), 
 * we'll walk through the list of "leaves" and check each one
 * for how many snps fall inside each. If it is larger than
 * the threshold, we'll replace it with it's geneBins for each
 * of the genes. 
 * 
 * If the user wishes to drill further yet, we'll check each
 * of the genes for the presence of 2 or more subgroups and
 * replace them with their intronic/exonic oriented bins.
 * 
 * Finally, if the user wishes to express by functionality, 
 * each bin will be evaluated for the number of function based
 * bins and the bin entry will be replaced by 2 or more
 * functional bins
 * 
 * At any time of expansion/expression, we'll just drop a 
 * bin if it's contents falls below the suggested minimum
 * threshold.
 * 
 * Once the bin list(s) are created, we'll simply construct
 * the snpIndex->bin lookup and destroy the intermediate bin
 * stuff. 
 * 
 * Locus Remap shows what the current index is and what it was. 
 * That is the index we actually store the index->bin/genotype
 * information for
 */
void BinManager::InitBins(
		const map<uint, Knowledge::GroupCollection*> &groups,
		const Knowledge::RegionCollection& regions,
		const vector<Knowledge::Locus*>& loci) {
	/**
	 * The following items will result in a bin once we reach the end of the 
	 * function. So, if we refine an item at one level, we should remove it from
	 * the original structure. These will be the index into the dataset-not into
	 * the vcf file (the locusRemap is used at the very last step only)
	 */
	/** Group -> snpIdx*/

	vector<Knowledge::Locus*>::const_iterator l_itr = loci.begin();
	vector<Knowledge::Locus*>::const_iterator l_end = loci.end();

	_total_variants = loci.size();

	while(l_itr != l_end){
		Knowledge::Locus& l = **l_itr;
		if ((1.0 - l.majorAlleleFreq()) < mafCutoff) {

			_rare_variants.insert(&l);

			// First, find all of the regions that contain this locus
			RegionCollection::const_region_iterator r_itr =
					regions.positionBegin(l.getChrom(), l.getPos());
			RegionCollection::const_region_iterator r_end =
					regions.positionEnd(l.getChrom(), l.getPos());

			Bin* curr_bin;
			if (r_itr == r_end){
				//Add to intergenic
				pair<short, int> key = make_pair(l.getChrom(), l.getPos() / IntergenicBinWidth);
				map<pair<short, int>, Bin*>::const_iterator i_bin =
						_intergenic_bins.find(key);

				if (i_bin == _intergenic_bins.end()){
					// OK, this bin is nonexistent
					curr_bin = new Bin(key.first, key.second);
					_intergenic_bins.insert(make_pair(key, curr_bin));
					_bin_list.insert(curr_bin);
				}else{
					curr_bin = (*i_bin).second;
				}
				curr_bin->addLocus(&l);
				_locus_bins[&l].insert(curr_bin);
			}
			while (r_itr != r_end){
				// For each region, get all of the groups it is in
				Region::const_group_iterator g_itr = (*r_itr)->groupBegin();
				Region::const_group_iterator g_end = (*r_itr)->groupEnd();

				// If not in any groups, add to region bins
				if(g_itr == g_end){
					int id = (*r_itr)->getID();
					map<int, Bin*>::const_iterator rm_itr = _region_bins.find(id);
					if (rm_itr == _region_bins.end()){
						curr_bin = new Bin(*r_itr);
						_bin_list.insert(curr_bin);
						_region_bins.insert(make_pair(id, curr_bin));
					}else{
						curr_bin = (*rm_itr).second;
					}
					curr_bin->addLocus(&l);
					_locus_bins[&l].insert(curr_bin);
				}

				// add to all group bins that it is a member of
				while(g_itr != g_end){
					int id = (*g_itr)->getID();
					map<int, Bin*>::const_iterator gm_itr = _group_bins.find(id);
					if (gm_itr == _group_bins.end()){
						curr_bin = new Bin(*g_itr);
						_bin_list.insert(curr_bin);
						_group_bins.insert(make_pair(id,curr_bin));
					}else{
						curr_bin = (*gm_itr).second;
					}
					curr_bin->addLocus(&l);
					_locus_bins[&l].insert(curr_bin);
					++g_itr;
				}
				++r_itr;
			}
		}
		++l_itr;
	}

	// At this point, we have all of the top level bins constructed and
	// stored in the variable _bin_list.  We should now collapse the
	// bins according to the preferences given
	collapseBins();
}

void BinManager::printBins(std::ostream& os, Knowledge::Locus* l,
		const string& sep) const{
	map<Knowledge::Locus*, set<Bin*> >::const_iterator m_itr = _locus_bins.find(l);

	if(m_itr != _locus_bins.end()){
		set<Bin*>::const_iterator s_itr = (*m_itr).second.begin();
		set<Bin*>::const_iterator s_end = (*m_itr).second.end();
		if (s_itr != s_end){
			os << (*s_itr)->getName();
			while(++s_itr != s_end){
				os << sep << (*s_itr)->getName();
			}
		}
	}

}

void BinManager::collapseBins(){
	return;
}

} // namespace BioBin;





