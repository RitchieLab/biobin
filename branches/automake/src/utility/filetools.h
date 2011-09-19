/* 
 * File:   filetools.h
 * Author: torstees
 *
 * Created on March 28, 2011, 4:33 PM
 */

#ifndef FILETOOLS_H
#define	FILETOOLS_H
#include <fstream>
#include <string>

namespace Utility {


std::string ExtractBaseFilename(const char *filename);

bool FileExists(const char *filename);




inline
bool FileExists(const char *filename) {
	std::ifstream testFile(filename);
	return testFile.is_open();
}
}

#endif	/* FILETOOLS_H */

