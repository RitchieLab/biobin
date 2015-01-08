/*
 * ICompressdFile.cpp
 *
 *  Created on: Apr 30, 2014
 *      Author: jrw32
 */

#include "ICompressedFile.h"

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

namespace BioBin {
namespace Utility {

ICompressedFile::ICompressedFile() : std::istream(0) {}

ICompressedFile::ICompressedFile(const char* fn, std::ios_base::openmode mode)
	: std::istream(0) {
	this->open(fn, mode);
}

void ICompressedFile::open(const char* fn, std::ios_base::openmode mode){
	int extPos = std::string(fn).find_last_of('.');
	std::string ext = std::string(fn).substr(extPos+1);

	bool isgz = (boost::iequals(ext, "gz") || boost::iequals(ext, "z"));
	bool isbz = boost::iequals(ext, "bz");

	_base_f.open(fn, mode | ((isgz || isbz) ? std::ios_base::binary : std::ios_base::in));

	if(_base_f.rdstate() != std::ios_base::failbit && (isgz || isbz)){

		if(isgz){
			_infile.push(boost::iostreams::gzip_decompressor());
		} else if(isbz){
			_infile.push(boost::iostreams::bzip2_decompressor());
		}

		_infile.push(_base_f);

		this->rdbuf(_infile.rdbuf());
	} else {
		this->rdbuf(_base_f.rdbuf());
		this->setstate(_base_f.rdstate());
	}

}


void ICompressedFile::close(){
	_infile.reset();
}


}

}
