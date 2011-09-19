/* 
 * File:   types.h
 * Author: torstees
 *
 * Created on March 28, 2011, 1:47 PM
 */

#ifndef TYPES_H
#define	TYPES_H

#include <sys/types.h>

#include <vector>
#include <string>
#include <set>

namespace Utility {
/**
 * Ideally, this is a set of indexes into the local snp data
 */
typedef std::set<uint> IdCollection;
typedef std::vector<std::string> StringArray;


#define MAX_LINE_LENGTH 600000

}
#endif	/* TYPES_H */

