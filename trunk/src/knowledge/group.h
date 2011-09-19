/* 
 * File:   group.h
 * Author: torstees
 *
 * Since we are just aggregating group and region children, I'm using this strictly
 * as a storage mechanism. The group manager will do all the real work
 *
 * Created on March 11, 2011, 9:03 AM
 */

#ifndef GROUP_H
#define	GROUP_H


#include <vector>
#include "def.h"

namespace Knowledge {

struct Group {
public:
	Group(uint id = 0) : id(id), name(""), desc("") { }
	Group(uint id, std::string name, std::string desc="") : id(id), name(name), desc(desc) {}
	Group(const Group& orig) : id(orig.id), name(orig.name), desc(orig.desc), regions(orig.regions), groups(orig.groups) { }
	virtual ~Group() {}

	uint id;
	std::string name;
	std::string desc;
	Utility::IdCollection regions;
	Utility::IdCollection groups;
};

}

#endif	/* GROUP_H */

