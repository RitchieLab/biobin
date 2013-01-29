#ifndef UTILITY_STRINGS_H
#define UTILITY_STRINGS_H

#include <string>
#include <vector>
#include "types.h"
#include <sstream>
#include <iostream>

namespace LDUtility {

#ifdef WIN32 
#define DIR_SLASH "\\"
#else
#define DIR_SLASH "/"
#endif

using namespace std;

/**
 * @brief The following are a set of string manipulation functions that might be useful to lots of folks
 */

	/**
	 * @brief Extracts the filename (right of the rightmost "/" and left of the leftmost "."
	 */
	string ExtractBaseFilename(const char *filename);

	/**
	 * @brief attempts to break the filename into 3 components
	 */
	void SplitIntoComponents(const char *filename, string &path, string &basename, string &ext);


	/**
	 * @brief Extracts the filename leaving off any path information)
	 */
	string ExtractFilename(const char *maybeHasFullPath);


	/**
	 * @brief Extract tokens from the string, origString, based on delimiters, sep and puts them into tokens
	 */
	uint TokenizeString(const char *origString, vector<string>& tokens, const char *sep);
	uint ExtractRsNumber(const char *rsid);
	/**
	 * @brief Strip whitespace from end of a string
	 */
	string StripTrailingWhitespace(const char *word);
	
	string ToString(int val);
	string ToString(uint val);
	string ToString(double val, int prec=0);

	uint CountColumns(const char *line);
	
	string ExtractExtension(const char *filename);

	string StripExtension(const char *originalFilename);

	string EscapeSpaces(const char *filename);

	string StripQuotes(const char *filename);

	string ParseFilename(istream &s, const char *desc = "filename");

	string FileToString(const char *filename, const char *sep = ",");
	vector<string> FileToVector(const char *filename);
	/**
	 * @Brief Converts chromosome IDs into numbers (X, Y, MT become 23, 24, 25)
	 */
	int ChromToInt(const char *chrom);
	string IntToChrom(int chrom);
	/**
	 * @brief Quick check for a valid filename
	 */
	bool FileExists(const char *filename);

	string ThousandsFormat(double d);

	/**
	 * @Brief Convert a string to all caps
	 */
	string ToUpper(const char *mixedcase);


	template<typename T>
	string ToString(T& ar, const char *sep=",") {
		stringstream s;
		typename T::iterator itr = ar.begin();
		typename T::iterator end = ar.end();
		uint words = 0;
		while (itr != end) {
			if (words++ > 0)
				s<<sep;
			s<<*itr++;
		}
		return s.str();
	}


}


#endif //UTILITY_STRINGS_H
