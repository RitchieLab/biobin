/* 
 * File:   filetools.cpp
 * Author: torstees
 * 
 * Created on March 28, 2011, 4:33 PM
 */

#include "filetools.h"
#include <string>
namespace Utility {
const char *seperators = "\\/";
/**
 * @brief Extracts the filename (right of the rightmost "/" and left of the leftmost "."
 */
std::string ExtractBaseFilename(const char *originalFilename) {
	std::string filename(originalFilename);

	size_t lastSlash = filename.find_last_of(seperators);
	if (lastSlash != std::string::npos)
		filename.erase(0, lastSlash + 1);

	size_t lastDot = filename.find_last_of(".");
	if (lastDot != std::string::npos)
		filename.erase(lastDot);

	return filename;
}



}
