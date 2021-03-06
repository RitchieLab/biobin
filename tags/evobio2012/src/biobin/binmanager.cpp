/* 
 * File:   binmanager.cpp
 * Author: torstees
 * 
 * Created on July 21, 2011, 3:54 PM
 */

#include "binmanager.h"

#include "PopulationManager.h"

#include "knowledge/RegionCollection.h"
#include "knowledge/Region.h"

using Knowledge::RegionCollection;
using Knowledge::Region;
using Knowledge::Group;
using std::make_pair;

namespace BioBin {

uint BinManager::IntergenicBinWidth = 50;
uint BinManager::BinTraverseThreshold = 50;
uint BinManager::MinBinSize = 1;
bool BinManager::ExpandByExons = false;
bool BinManager::ExpandByFunction = false;
float BinManager::mafCutoff = 0.05;
uint BinManager::maxSnpCount = 200;

BinManager::BinManager(const PopulationManager& pop_mgr) : _total_variants(0),
		_pop_mgr(pop_mgr){}

BinManager::~BinManager() {
	set<Bin*>::const_iterator itr = _bin_list.begin();
	set<Bin*>::const_iterator end = _bin_list.end();
	while(itr != end){
		delete *itr;
		++itr;
	}
}

void BinManager::InitBins(
		const map<uint, Knowledge::GroupCollection*> &groups,
		const Knowledge::RegionCollection& regions,
		const vector<Knowledge::Locus*>& loci) {

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
				pair<short, int> key = make_pair(l.getChrom(), l.getPos() / (IntergenicBinWidth*1000));
				map<pair<short, int>, Bin*>::const_iterator i_bin =
						_intergenic_bins.find(key);

				if (i_bin == _intergenic_bins.end()){
					// OK, this bin is nonexistent
					curr_bin = new Bin(_pop_mgr, key.first, key.second);
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
					curr_bin = addRegionBin(*r_itr);
					curr_bin->addLocus(&l);
					_locus_bins[&l].insert(curr_bin);
				}

				// add to all group bins that it is a member of
				while(g_itr != g_end){
					int id = (*g_itr)->getID();
					map<int, Bin*>::const_iterator gm_itr = _group_bins.find(id);
					if (gm_itr == _group_bins.end()){
						curr_bin = new Bin(_pop_mgr, *g_itr);
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
	// First things first, let's expand the groups that are too big
	set<Bin*>::iterator b_itr = _bin_list.begin();
	Bin::const_locus_iterator v_itr;
	Bin::const_locus_iterator v_end;
	while(b_itr != _bin_list.end() && (*b_itr)->isGroup()){
		if((uint)(*b_itr)->getSize() > BinTraverseThreshold){

			// First, add all of the appropriate child bins
			Group* curr_group = (*b_itr)->getGroup();
			Group::const_region_iterator r_itr = curr_group->regionBegin();
			Group::const_region_iterator r_end = curr_group->regionEnd();



			Bin* curr_bin;
			v_end = (*b_itr)->variantEnd();
			while(r_itr != r_end){

				// For all rare loci in the current region, associate them
				// with this bin
				v_itr = (*b_itr)->variantBegin();

				while(v_itr != v_end){
					if ((*r_itr)->containsLocus(**v_itr)){
						curr_bin = addRegionBin(*r_itr);
						curr_bin->addLocus(*v_itr);
						_locus_bins[*v_itr].insert(curr_bin);
					}
					++v_itr;
				}
				++r_itr;
			}

			// Now, erase the bin
			eraseBin(b_itr);
		}else{
			++b_itr;
		}
	}

	//OK, now we go through and clean up all of the bins that are too small!
	b_itr = _bin_list.begin();
	while(b_itr != _bin_list.end()){
		if((uint)(*b_itr)->getSize() < MinBinSize){
			eraseBin(b_itr);
		}else{
			// Otherwise, just move along, nothing to see here
			++b_itr;
		}
	}
	return;
}

void BinManager::eraseBin(set<Bin*>::iterator& b_itr){

	// Remove the bin from the shortcut data structures
	if((*b_itr)->isGroup()){
		_group_bins.erase((*b_itr)->getID());
	}else if((*b_itr)->isIntergenic()){
		_intergenic_bins.erase(
				make_pair((*b_itr)->getChrom(), (*b_itr)->getID()));
	}else{
		_region_bins.erase((*b_itr)->getID());
	}

	// remove the bin from the locus->set<Bin*>
	Bin::const_locus_iterator v_itr = (*b_itr)->variantBegin();
	Bin::const_locus_iterator v_end = (*b_itr)->variantEnd();
	while(v_itr != v_end){
		_locus_bins[*v_itr].erase(*b_itr);
		++v_itr;
	}

	// Delete the bin itself
	delete *b_itr;

	// And remove it from the master list and move to the next one
	_bin_list.erase(b_itr++);
}

Bin* BinManager::addRegionBin(Region* reg){
	Bin* curr_bin;
	int id = reg->getID();
	map<int, Bin*>::const_iterator rm_itr = _region_bins.find(id);
	if (rm_itr == _region_bins.end()){
		curr_bin = new Bin(_pop_mgr, reg);
		_bin_list.insert(curr_bin);
		_region_bins.insert(make_pair(id, curr_bin));
	}else{
		curr_bin = (*rm_itr).second;
	}
	return curr_bin;
}

} // namespace BioBin;





