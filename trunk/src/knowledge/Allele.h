/*
 * Allele.h
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_ALLELE_H
#define KNOWLEDGE_ALLELE_H

#include <string>

using std::string;

namespace Knowledge{

class Allele{
public:
	Allele(const string& data, float freq) : _data(data), _freq(freq){}
	~Allele() {}

	bool operator<(const Allele&) const;
	bool operator>(const Allele&) const;
	bool operator==(const Allele&) const;

	float getFreq() const { return _freq;}
	// Will we need to access the allele?  we'll see
	const string& getData() const { return _data;}

private:

	string _data;
	float _freq;
};

}

#endif /* ALLELE_H_ */
