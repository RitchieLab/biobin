/*
 * TestFactory.h
 *
 *  Created on: Feb 12, 2015
 *      Author: jrw32
 */

#include <map>
#include <string>

#ifndef TESTFACTORY_H_
#define TESTFACTORY_H_

namespace BioBin {

namespace Test {

class Test;

typedef Test* (createFunc)();

class TestFactory {
public:

	typedef std::map<const std::string, createFunc*>::const_iterator const_iterator;

private:
	TestFactory(){}
	TestFactory(const TestFactory&);
	TestFactory& operator=(const TestFactory&);

public:
	const std::string& RegisterTest(const std::string& key, createFunc* ptr);
	Test* Create(const std::string& key);

	const_iterator begin() const{return creation_map.begin();}
	const_iterator end() const{return creation_map.end();}

	static TestFactory& getFactory(){static TestFactory f; return f;}


private:
	std::map<const std::string, createFunc*> creation_map;

};

}

}

#endif /* TESTFACTORY_H_ */
