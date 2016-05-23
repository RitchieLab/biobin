/*
 * Phenotype.h
 *
 *  Created on: Feb 17, 2015
 *      Author: jrw32
 */

#ifndef BIOBIN_UTILITY_PHENOTYPE_H
#define BIOBIN_UTILITY_PHENOTYPE_H

#include <boost/dynamic_bitset.hpp>

namespace BioBin {
namespace Utility {

class Phenotype {
public:
	typedef std::pair<boost::dynamic_bitset<>, boost::dynamic_bitset<> > bitset_pair;

	Phenotype(unsigned int idx, const bitset_pair& status) :
		_idx(idx), _status(&status) {
	}
	const bitset_pair& getStatus() const {
		return *_status;
	}
	unsigned int getIndex() const {
		return _idx;
	}

private:
	unsigned int _idx;
	//! bistset pair.  Note: _status.first == control, _status.second==case.
	const bitset_pair* _status;
};

}
}

#endif /* PHENOTYPE_H_ */
