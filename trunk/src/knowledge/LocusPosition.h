/*
 * LocusPostion.h
 *
 *  Created on: Jan 13, 2015
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_LOCUSPOSTION_H
#define KNOWLEDGE_LOCUSPOSTION_H

namespace Knowledge {

class LocusPosition {
public:
	LocusPosition(unsigned short chr, unsigned int pos){
		_data._chr = chr & UNKNOWN_CHR;
		_data._pos = pos;
	}
	~LocusPosition() {}

	unsigned short getChrom() const { return _data._chr == UNKNOWN_CHR ? UNKNOWN_CHR_RETURN : _data._chr;}
	unsigned int getPos() const {return _data._pos;}

	// default copy/assignment ctors OK here!

private:
	struct chrpos{
		unsigned int _chr : 5;
		unsigned int _pos : 27;
	} _data;

	static const unsigned short UNKNOWN_CHR=31;
	static const unsigned short UNKNOWN_CHR_RETURN=static_cast<unsigned short>(-1);
};

}

#endif /* LOCUSPOSTION_H_ */
