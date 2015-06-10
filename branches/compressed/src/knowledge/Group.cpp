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

}



