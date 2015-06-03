/*
 * ICompressdFile.h
 *
 *  Created on: Apr 30, 2014
 *      Author: jrw32
 */

#ifndef UTILITY_ICOMPRESSEDFILE_H
#define UTILITY_ICOMPRESSEDFILE_H

#include <istream>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>

namespace BioBin {
namespace Utility {

/*!
 * \brief A class to read in compressed files automatically based on extension
 * This class, which has the same interface as an ifstream (so can be used
 * interchangeably) will automatically decompress
 */
class ICompressedFile : public std::istream {
public:
	ICompressedFile();
	explicit ICompressedFile(const char* fn, std::ios_base::openmode mode = std::ios_base::in);
	virtual ~ICompressedFile() {}

	void open(const char* fn, std::ios_base::openmode mode = std::ios_base::in);
	void close();

private:
	void openCompressed(const char* fn);

	boost::iostreams::filtering_istream _infile;
	std::ifstream _base_f;


};

}

}

#endif /* ICOMPRESSDFILE_H_ */
