/*
 * TestFactory.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include "TestFactory.h"

using std::map;
using std::string;

namespace BioBin {

namespace Test {

const string& TestFactory::RegisterTest(const string& key, createFunc* ptr){
	creation_map[key] = ptr;
	return key;
}

Test* TestFactory::Create(const string& key){
	map<string, createFunc*>::const_iterator it=creation_map.find(key);
	if(it != creation_map.end()){
		return (*it).second();
	}else{
		return NULL;
	}
}

}

}
