/*
 * Information.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: jrw32
 */

#include "Information.h"

namespace Knowledge{

const Information::snp_role Information::EXON("exon");
const Information::snp_role Information::INTRON("intron");
const Information::snp_role Information::REGULATORY("reg");

map<string, int> Information::snp_role::s_val_map;
set<Information::snp_role*, Information::snp_role::Ptr_Less> Information::snp_role::s_enums;
int Information::snp_role::s_num_vals = 0;

}



