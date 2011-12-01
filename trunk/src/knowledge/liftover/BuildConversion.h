/*
 * BuildConversion.h
 *
 *  Created on: Nov 30, 2011
 *      Author: jrw32
 */

#ifndef KNOWLEDGE_LIFTOVER_BUILDCONVERSION_H
#define KNOWLEDGE_LIFTOVER_BUILDCONVERSION_H

#include <string>

using std::string;

namespace Knowledge{
namespace Liftover{

class BuildConversion {
public:
	BuildConversion(const string& chrom, int start, int stop);

	void align();
	bool operator<(const BuildConversion& other) const;

	void setStarts(int local, int remote){
		_lStart = local;
		_rStart = remote;
	}
	void setStops(int local, int remote){
		_lStop = local;
		_rStop = remote;
	}
	void setChroms(const string& local, const string& remote){
		_lChrom = local;
		_rChrom = remote;
	}
	void setLocalStart(int local){
		_lStart = local;
	}
	void setLocalStop(int local){
		_lStop = local;
	}

	int getLocalStart() const { return _lStart;}
	int getLocalStop() const {return _lStop;}
	int getRemoteStart() const {return _rStart;}
	int getRemoteStop() const {return _rStop;}
	const string& getLocalChrom() const {return _lChrom;}
	const string& getRemoteChrom() const {return _rChrom;}

private:
	long _score;					///< The conversion's score
	string _lChrom;		///< Local Chromosome
	int _lStart;					///< Local Start
	int _lStop;					///< Local Stop
	string _rChrom;		///< Remote Chromosome
	int _rStart;					///< Remote Start
	int _rStop;					///< Remote Stop
};

}
}

#endif /* BUILDCONVERSION_H_ */
