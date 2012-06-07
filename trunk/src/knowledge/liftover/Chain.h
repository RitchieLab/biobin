/*
 * Chain.h
 *
 *  Created on: May 18, 2012
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_LIFTOVER_CHAIN_H
#define KNOWLEDGE_LIFTOVER_CHAIN_H

#include <set>
#include <utility>

#include "Segment.h"

using std::pair;
using std::make_pair;
using std::set;

namespace Knowledge{

namespace Liftover{

class Chain {
public:

	enum error_codes{
		NOT_INTERSECTING,
		DELETED,
		PARTIALLY_DELETED
	};

public:
	Chain(int id, long score, int o_s, int o_e, short n_c, bool fwd) :
		_id(id), _score(score), _old_start(o_s), _old_end(o_e),
		_new_chrom(n_c), _is_fwd(fwd) {};


	bool overlaps(int start, int end) const {
		return _old_start <= end && _old_end >= start;
	}
	bool operator<(const Chain& other) const {
		return (_score == other._score) ? _id < other._id : _score > other._score;
	}

	short getNewChrom() const{
		return _new_chrom;
	}

	int getID() const{
		return _id;
	}

	pair<int, int> convertRegion(int start, int end, float minMappingFrac=0.95) const;

	void addSegment(int old_s, int old_e, int new_s);

private:
	// No copying or assignment, please (though it would probably be OK)
	Chain(const Chain& other);
	Chain& operator=(const Chain& other);


	int _id;
	long _score;
	int _old_start;
	int _old_end;
	short _new_chrom;
	bool _is_fwd;

	set<Segment> _data;
};

}
}

namespace std{

template<>
struct less<Knowledge::Liftover::Chain*> {

	bool operator()(const Knowledge::Liftover::Chain* x, const Knowledge::Liftover::Chain* y) const{
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};

}

#endif /* CHAIN_H_ */
