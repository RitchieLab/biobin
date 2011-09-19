/* 
 * File:   strings.h
 * Author: torstees
 *
 * Created on March 8, 2011, 1:56 PM
 */

#ifndef STRINGS_H
#define	STRINGS_H

#include <string>
#include <fstream>
#include <vector>
#include <set>
#include "types.h"
#include <boost/tokenizer.hpp>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "filetools.h"


namespace Utility {

std::string UpperCase(const char *value);

std::string ENV(const char *var, const char *def);


/**
 * Loads contents of file into the string
 * @param filename
 * @return
 */
std::string LoadContents(const char *filename);


/**
 * Splits string on one or more separating characters-when strict is true, empty strings are returned (i.e. a,,b is ['a','','b'}
 * @param s The string being split
 * @param sepChars The seperation characters
 * @param strict T/F indicating whether or not to allow empty strings (default, false, returns no empty strings)
 * @return
 */
StringArray Split(const char *source, const char *delimiters = " \t\n", bool keep_empty_tokens = false);


/**
 * Stringsplut differs from Split because it interpresents delimiter as a word, so you can
 * split on words, rather than one or more characters
 * @param source
 * @param delimiters
 * @param keep_empty_tokens
 * @return
 */
StringArray StringSplit(const char *source, const char *delimiter = ",", bool keep_empty_tokens = false);


/**
 * This might not be the most efficient approach, since it uses stringstreams to
 * perform the conversion from type to type. 
 * @param source
 * @param delimiters
 * @param keep_empty_tokens
 * @return
 */
template<typename T>
std::vector<T> Split(const char *source, const char *delimiters = " \t\n", bool keep_empty_tokens = false);

/*
template <typename T>
std::string Join(const T& data, const char *delimiter = "\t");
*/

template <typename T>
std::string Join(const T& data, const char *delimeter="\t", uint start = 0, uint stop = -1);

template <typename T>
std::string MapJoin(const T& data, const char *delimeter="\t", const char *kep_sep=":");

template <typename T>
std::set<T> ToSet(const char *s, const char *del = " \t\n");

template<typename T>
std::string ToString(T& ar, const char *sep);

template<typename T>
std::string ToString(T val);

/**
 * Joins the string array and returns the string (uses delimiter to join)
 * @param strings
 * @param delimiter
 * @return
 */
std::string Join(const StringArray& strings, const char *delimiter = "\t");


std::string GetPartialString(const char *basestring, uint start, uint stop=-1, const char *del = "\t");


std::string StripTrailingWhitespace(const char *word);

template <typename T>
inline std::string ToString(T val) {
	std::stringstream ss;
	ss<<val;
	return ss.str();
}

inline std::string GetPartialString(const char *basestring, uint start, uint stop, const char *del) {
	StringArray strings = Split(basestring, del);
	return Join(strings, del, start, stop);
}

template <typename T>
inline std::set<T> ToSet(const char *s, const char *del) {
	std::set<T> output;
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(del);
	std::string source(s);

	tokenizer tokens(source, sep);
	for (tokenizer::iterator itr=tokens.begin(); itr!=tokens.end(); itr++) {
		T val;
		std::stringstream ss(*itr);
		ss>>val;
		output.insert(val);
	}
	return output;
}
template <typename T>
inline std::string MapJoin(const T& data, const char *del="\t", const char *key_sep=":") {
	typename T::const_iterator itr = data.begin();
	typename T::const_iterator end = data.end();

	uint count = 0;
	std::stringstream ss;
	while (itr != end) {
			if (count++)
				ss<<del;
			ss<<itr->first<<key_sep<<itr->second;
		itr++;
	}
	return ss.str();
}
template <typename T>
inline std::string Join(const T& data, const char *del, uint start, uint stop) {
	typename T::const_iterator itr = data.begin();
	typename T::const_iterator end = data.end();

	std::stringstream ss;
	uint count = 0;
	uint index = 0;

	while (itr != end) {
		if (index >= start && (index < stop || index == (uint)-1) ) {
			if (count++)
				ss<<del;
			ss<<*itr;
		}
		itr++;
		index++;
	}
	return ss.str();
}

template <typename T>
std::vector<T> Split(const char *s, const char *del, bool keep_empty_tokens) {
	std::vector<T> output;
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(del);
	boost::char_separator<char> sepstrict(del, "", boost::keep_empty_tokens);
	std::string source(s);

	tokenizer tokens(source, sep);
	if (keep_empty_tokens)
		tokens = tokenizer(source, sepstrict);
	for (tokenizer::iterator itr=tokens.begin(); itr!=tokens.end(); itr++) {
		T val;
		std::stringstream ss(*itr);
		ss>>val;
		output.push_back(val);
	}
	return output;
}

/*

template <typename T>
inline std::string Join(const T& data, const char *delimiter) {
	typename T::const_iterator itr = data.begin();
	typename T::const_iterator end = data.end();
	std::stringstream ss;
	int count = 0;
	while (itr != end) {
		if (count++>0)
			ss<<delimiter;
		ss<<*itr++;
	}
	return ss.str();
}
 * */

inline
std::string ENV(const char *var, const char *def) {
	char *variable = getenv(var);
	if (variable)
		return variable;
	return def;
}

inline
std::string Join(const StringArray& strings, const char *delimiter) {
	return Join<StringArray>(strings, delimiter, 0, -1);
}

inline
StringArray StringSplit(const char *s, const char *del, bool keep_empty_tokens) {
	StringArray strings;
	std::string base = s;
	uint delLen = strlen(del);
	uint start = 0;
	size_t cur = base.find(del);

	while (cur != std::string::npos) {
		std::string token = base.substr(start, cur-start);
		if (keep_empty_tokens || token.length() > 0)
			strings.push_back(token);
		start = cur + delLen;
		cur = base.find(del, start);
	}

	std::string token = base.substr(start, cur);
	if (keep_empty_tokens || token.length() > 0)
		strings.push_back(token);
	return strings;
}

inline
StringArray Split(const char *s, const char *del, bool keep_empty_tokens) {
	StringArray strings;
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(del);
	boost::char_separator<char> sepstrict(del, "", boost::keep_empty_tokens);
	std::string source(s);
	
	tokenizer tokens(source, sep);
	if (keep_empty_tokens)
		tokens = tokenizer(source, sepstrict);
	for (tokenizer::iterator itr=tokens.begin(); itr!=tokens.end(); itr++)
		strings.push_back(*itr);
	return strings;
}

inline
std::string LoadContents(const char *filename) {
	std::ifstream file(filename, std::ios::binary);

	std::streampos length = file.tellg();
	file.seekg(0,std::ios::end);
	length = file.tellg() - length;
	
	std::string contents;

	if ((long)length > 0 ) {
		file.seekg(0, std::ios::beg);

		char *buffer;
		buffer = (char*) malloc((long)length+1);
		file.read(buffer, (long)length);
		buffer[(long)length] = 0;
		file.close();
		contents = buffer;
		free(buffer);

		//Handle DOS carriage returns
		char badChar[2];
		badChar[0] = char(13);
		badChar[1] = 0;
		boost::replace_all(contents, std::string(badChar), "");
	}
	return contents;
}


inline
std::string UpperCase(const char *lower) {
	std::string s(lower);
	uint count = s.length();
	for (uint i=0; i<count; i++) 
		s[i] = toupper(s[i]);
	return s;
}

template<typename T>
inline std::string ToString(T& ar, const char *sep=",") {
	std::stringstream s;
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


// Quick way to remove whitespace from strings (windows /r for the most part from end of lines)
// http://www.cplusplus.com/reference/string/string/find_last_not_of/
inline
std::string StripTrailingWhitespace(const char *word){ 
	std::string whitespace("\t\f\v\n\r ");
	std::string str = word;
	size_t found=str.find_last_not_of(whitespace);
	if (found!=std::string::npos)
		return str.erase(found+1);
	else 
		return "";
}

}
#endif	/* STRINGS_H */

