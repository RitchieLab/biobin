/*
 * Segment.h
 *
 *  Created on: May 18, 2012
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_LIFTOVER_SEGMENT_H
#define KNOWLEDGE_LIFTOVER_SEGMENT_H

namespace Knowledge{

namespace Liftover {

class Segment {
public:
	Segment(int old_s, int old_e, int new_s) :
		_old_start(old_s), _old_end(old_e), _new_start(new_s) {};

	/**
	 * \brief Implicit conversion from int.
	 * This constructor provides an implicit conversion from integer -> Segment
	 * that is necessary when finding all segments whose start value is compared
	 * to a random integer.
	 *
	 * \param The starting value to use as a searching key
	 */
	Segment(int start_key=0) :
		_old_start(start_key), _old_end(start_key), _new_start(start_key) {};

	// NOTE: compiler-generated copy/assignment/destructor are FINE!

	int size() const {
		return _old_end - _old_start;
	}

	int getNewStart() const {
		return _new_start;
	}

	int getOldStart() const {
		return _old_start;
	}

	int getOldEnd() const {
		return _old_end;
	}

	bool operator<(const Segment& other) const {
		return _old_start < other._old_start;
	}



private:
	int _old_start;
	int _old_end;
	int _new_start;
};

}
}


#endif /* SEGMENT_H_ */
