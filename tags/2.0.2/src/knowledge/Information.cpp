/*
 * Information.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: jrw32
 */

#include "Information.h"

#include <sstream>

using std::stringstream;
using std::ostream;
using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;
using boost::unordered_map;

namespace Knowledge{
map<string, unsigned long> Information::snp_role::s_val_map;
set<const Information::snp_role*, Information::snp_role::Ptr_Less> Information::snp_role::s_enums;
int Information::snp_role::s_num_vals = 0;

const Information::snp_role Information::EXON("exon");
const Information::snp_role Information::INTRON("intron");
const Information::snp_role Information::REGULATORY("reg");
const Information::snp_role Information::OTHER("other");

vector<string> Information::c_source_names;
vector<string> Information::c_role_files;
vector<string> Information::c_source_exclude;
set<unsigned int> Information::_s_source_ids;

string Information::getSourceList(){
	string ret_val = "";
	//if(c_source_names.size() > 0){
		const set<unsigned int>& id_set = getSourceIds();
		if(id_set.size() > 0){
			stringstream id_str;
			id_str << "(";
			set<unsigned int>::const_iterator itr = id_set.begin();
			id_str << *itr;
			while(++itr != id_set.end()){
				id_str << "," << *itr;
			}
			id_str << ")";

			ret_val = id_str.str();
		}
	//}

	return ret_val;
}

}



