/*
 * Group.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: jrw32
 */

#include "Group.h"

using std::string;
using std::set;

namespace Knowledge{

Group::Group(uint id, const string& name, const string& desc) :
		_id(id), _name(name), _description(desc) {}


Group::~Group(){
	// If I'm deleting this group, make sure to get rid of dangling pointers
	for(set<Region*>::const_iterator ri=_regions.begin(); ri != _regions.end(); ri++){
		(*ri)->removeGroup(*this);
	}
	for(set<Group*>::const_iterator gi=_parents.begin(); gi!=_parents.end(); gi++){
		(*gi)->_children.erase(this);
	}
	for(set<Group*>::const_iterator gi=_children.begin(); gi!=_children.end(); gi++){
		(*gi)->_parents.erase(this);
	}

}

}



