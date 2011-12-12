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
#include <stdlib.h>

using std::string;

using std::ostream;

namespace Knowledge{

class Allele{
public:
	Allele(const string& data, float freq, uint pos) : _data(data), _freq(freq), _pos(pos){}
	~Allele() {}

	bool operator<(const Allele&) const;
	bool operator>(const Allele&) const;
	bool operator==(const Allele&) const;

	float getFreq() const { return _freq;}
	// Will we need to access the allele?  we'll see
	const string& getData() const { return _data;}
	// Returns the position of the allele (0 is reference typically)
	uint getPos() const {return _pos;}

	void print(ostream& o, const string& sep=":") const;

private:

	string _data;
	float _freq;
	uint _pos;
};

}

ostream& operator<<(ostream& o, const Knowledge::Allele& l);

#endif /* ALLELE_H_ */
