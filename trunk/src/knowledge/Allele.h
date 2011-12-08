/*
 * Allele.h
 *
 *  Created on: Nov 15, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_ALLELE_H
#define KNOWLEDGE_ALLELE_H

#include <string>
#include <ostream>

using std::string;

using std::ostream;

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

	void print(ostream& o, const string& sep=":") const;

private:

	string _data;
	float _freq;
};

}

ostream& operator<<(ostream& o, const Knowledge::Allele& l);

#endif /* ALLELE_H_ */
