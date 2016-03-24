/*
 * Wilcoxon.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: jrw32
 */

#include "Wilcoxon.h"

#include <vector>
#include <utility>
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#include <gsl/gsl_cdf.h>

using std::vector;
using std::pair;
using std::string;

namespace BioBin {

namespace Test {

string Wilcoxon::testname = Wilcoxon::doRegister("wilcoxon");

void Wilcoxon::init(){
	if(_pop_mgr_ptr->getNumCovars() > 0){
		std::cerr << "WARNING: The Wilcoxon test ignores covariates, which were given!" << std::endl;
	}
}

double Wilcoxon::runTest(const Bin& bin) const{

	// let's get the values for all non-missing samples
	const Utility::Phenotype::bitset_pair& status = _pheno_ptr->getStatus();

	boost::dynamic_bitset<> nonmiss = status.first | status.second;

	vector<std::pair<float, unsigned int> > data;
	data.reserve(status.first.size());

	// case/control status doesn't really matter here, as long as it's not missing!
	for(unsigned int i=0; i<status.first.size(); i++){
		if(nonmiss[i]){
			data.push_back(std::make_pair(_pop_mgr_ptr->getTotalIndivContrib(bin, i,  *_pheno_ptr), i));
		}
	}

	// now sort the data from smallest to largest
	std::sort(data.begin(), data.end());

	double W = 0;
	//iterate over the data, and add the rank to W if status.first[pos] is set
	// TODO: account for ties
	vector<unsigned int> numtie;
	unsigned int i=0;
	while(i<data.size()){
		unsigned int st=i;
		float currv = data[i].first;
		while(i<data.size() - 1 && data[i+1].first == currv){++i;}
		++i;
		// now, st will be where I started and i will be the next value
		if(st < i-1){
			numtie.push_back(i-st-1);
		}

		// If I'm here, then I need to take an average
		double avg = (i-st+1)/2.0 + st;

		for(unsigned int j=st; j<i; j++){
			if(status.first[data[j].second]){
				W += avg;
			}
		}
	}

	// now, construct the Z statistic
	unsigned int m=status.first.count();
	unsigned int n=status.second.count();

	double denom = m*n*(m+n+1)/12.0;
	double num = W - m*(m+n+1)/2.0;

	double adj = 0;
	for(unsigned int j=0; j<numtie.size(); j++){
		adj+= (numtie[j]-1)*numtie[j]*(numtie[j]+1);
	}
	if(adj != 0){
		denom -= m*n/(12.0*(m+n)*(m+n-1))*adj;
	}


	double Z = num / sqrt(denom);

	// now check the Z statistic, taking the area from the std normal
	double pval;
	if(Z > 0){
		pval = 2*gsl_cdf_ugaussian_Q(Z);
	} else {
		pval = 2*gsl_cdf_ugaussian_P(Z);
	}

	return pval;
}

}

}
