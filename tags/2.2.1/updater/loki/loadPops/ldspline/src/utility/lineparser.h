//
// C++ Interface: lineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESELINEPARSER_H
#define ESELINEPARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "types.h"
#include "strings.h"
#include <map>
#include <cstring>

namespace LDUtility {

using namespace std;


/**
 * @brief Baseclass for the ascii parser routines
 */
class AsciiParser {
public:
	AsciiParser(){}
	virtual ~AsciiParser() {}
	virtual bool ParseLine(const char *line, uint val)=0;
};


/**
  @brief Breaks the line into a vector and passes it to the functor
 */
class ArrayedTextParser : public AsciiParser {
public:
	ArrayedTextParser() {}
	virtual ~ArrayedTextParser() {} 
	bool ParseLine(const char *line, uint val);
	virtual bool Evaluate(StringArray& words, uint val)=0;
};



/**
@brief Parses ascii files and executes a functor for each line

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LineParser{
public:
    LineParser(char cmtChar = 0) : cmtChar(cmtChar) {}

    ~LineParser(){}

	/**
	 * @brief Opens a file and sends each line to the function pointer at fn
	 */
	uint Parse(const char* filename, AsciiParser *parser, bool keepComments = false);

	string comments;
protected:
	char cmtChar;
};

/**
  @brief Converts a file to an array of strings
 */
class FileToArray : public AsciiParser {
public:
	FileToArray();
	~FileToArray();
	
	bool ParseLine(const char *line, uint val=0);
	StringArray OpenFile(const char *filename,char cmtChar = 0);
	StringArray strings;
};



class FileToMap: public AsciiParser {
public:

	FileToMap(char comment = '#');
	virtual ~FileToMap();

	bool ParseLine(const char *line, uint val=0);


	/**
	 * @brief Opens file and parses it, returning the number of valie lines
	 * @param filename The file to be parsed
	 * @param commentChar If the first character matches this, the line is treated as a comment
	 * @param true indicates that all comments will be stored as a single large string
	 */
	uint Parse(const char *filename);
	/**
	 * @brief Returns the vector of lines associated with a given key. The key will remain as the first member of the line
	 * @param key The first word in the line
	 * @param lines a vector containing each string that was found
	 * @return t/f indicating that the key was found
	 * @note lines will not be changed if the return value is False. 
	 * It is important to note that the entire contents of the file will be loaded into memory.
	 */
	bool GetLines(const char *key, vector<string>& lines);

	double GetDouble(const char *key);
	char GetChar(const char *key);
	/**	
	 * @brief Sets value at key to value. 
	 * @note This assumes that there is only one value, so any prev. value(s) is removed
	 */
	void SetValue(const char *key, double value);

	virtual void AppendValue(const char *key, const char *value);

	int GetInteger(const char *key);
	/**	
	 * @brief Sets value at key to value. 
	 * @note This assumes that there is only one value, so any prev. value(s) is removed
	 */	
	void SetValue(const char *key, int value);

	void SetValue(const char *key, const string& value);
	
	string GetLine(const char *key);

	bool GetBoolean(const char *key);

	/**
	 * @Brief Returns a single string for all elements of the list
	 */
	string GetString(const char *key);

	vector<string> GetKeys() { return keys; }

	void Write(ostream& os);
	string comments;
protected:
	map<string, vector<string> > strings;
	map<string, string> varComments;
	vector<string> keys;
	char comment;


};

/**
 * @brief Converts file to a map, where each key MUST be found in the "keys" list or else determined an error
 */
class FileToMapExplicit: public FileToMap {
public:
	FileToMapExplicit(char comment = '#') : FileToMap(comment) { }
	bool ParseLine(const char *line, uint val=0);

	/**
	 * @brief Initialize a key, storing the default value and the comment
	 * @note Repeated calls work, except the comment that is recorded is the last
	 */
	void InitKey(const char *key, const char *defValue, const char *comment);

	/**
	 * @brief Parse the file
	 */
	uint Parse(const char *filename);

	void AppendValue(const char *key, const char *val);
protected:
	map<string, bool> defaultValue; 	///<Indicate that a given value has not been read from the file.
};


/**
 * @brief A simple parser object that returns the number of lines encountered. 
 * Eventually, this could be expanded to get other types of information like the number of columns in certain types of lines
 */
class BasicLineCounter: public AsciiParser {
public:
	BasicLineCounter() {}
	~BasicLineCounter() {}

	bool ParseLine(const char *line, uint val=0);
};


/**
@brief Sets values in the bitvector to value if they appear in the file
*/
class ExclusionList : public AsciiParser {
public:
	ExclusionList(BitSetType *active, bool value=false) : active(active), value(value) {}
	~ExclusionList() {}
	
	bool ParseLine(const char *line, uint val=0);
protected:
	BitSetType *active;
	bool value;

};


inline
bool ExclusionList::ParseLine(const char* line, uint val /*=0*/) {
	bool success = strlen(line) > 0;

	if (success)	{
		uint idx=atoi(line);	
		(*active)[idx-1]=value;
	}
	return success;
}



inline
bool ArrayedTextParser::ParseLine(const char *rawline, uint val) {
	StringArray line;
	size_t len=strlen(rawline);
	
	if (len > 0) {
		TokenizeString(rawline, line, " \t,\r");
		return Evaluate(line, val);
	}	
	return false;
}

inline
FileToMap::FileToMap(char comment) : comment(comment) { }

inline
FileToMap::~FileToMap() { }

inline
bool FileToMap::GetLines(const char *key, vector<string>& lines) {
	map<string, vector<string> >::iterator itr = strings.find(key);
	map<string, vector<string> >::iterator end = strings.end();

	bool isFound = false;
	if (itr != end) {
		isFound = true;
		lines=itr->second;
	}
	return isFound;
}

inline
string FileToMap::GetString(const char *key) {
	string value = "";
	vector<string> lines;
	
	if (GetLines(key,lines) && lines.size() > 0) {
		for (uint i=0; i<lines.size(); i++) {
			if (i > 0) 
				value += " ";	
			value += lines[i];
		}
	}
	return value;
}
inline
string FileToMap::GetLine(const char *key) {
	string value="";
	vector<string> lines;

	if (GetLines(key, lines) && lines.size()>0)
		value=lines[0];
	return value;
}

inline
bool FileToMap::GetBoolean(const char *key) {
	bool isOn=false;
	string v = GetLine(key);
	const char *value = v.c_str();
	isOn = strcmp(value, "ON")==0 || strcmp(value, "on")==0 || strcmp(value, "On")==0 || strcmp(value, "YES")==0 || strcmp(value, "Yes") ==0 || strcmp(value, "Yes")==0;
	return isOn;
}

inline
int FileToMap::GetInteger(const char *key) {
	return atoi(GetLine(key).c_str());
}

inline
char FileToMap::GetChar(const char *key) {
	string line = GetLine(key);
	if (line.length() > 0) 
		return line.c_str()[0];
	else
		return '\0';
}

inline
double FileToMap::GetDouble(const char *key) {
	return atof(GetLine(key).c_str());
}

inline
void FileToMap::AppendValue(const char *key, const char *val) {
	vector<string>& list = strings[key];
	list.push_back(val);
}

inline
void FileToMap::SetValue(const char *key, const string& val) {
	if (strings.find(key) != strings.end()) 
		strings[key] = vector<string>();
	vector<string>& list = strings[key];
	list.clear();
	list.push_back(val);
}

inline
void FileToMap::SetValue(const char *key, double value) {
	std::stringstream ss;
	ss<<value;
	string s=ss.str();
	SetValue(key, s);
}

inline
void FileToMap::SetValue(const char *key, int value) {
	std::stringstream ss;
	ss<<value;
	string s=ss.str();
	SetValue(key, s);
}

inline 
bool FileToMap::ParseLine(const char *line, uint val /* = 0*/) {
	stringstream source(line);

	bool valid = strlen(line) > 0 && line[0] != '\n';

	if (valid) {
		string key, word;

		source>>key;
	
		if (strings.find(key) == strings.end()) 
			keys.push_back(key);
		
		//We don't really care if contents is empty, 
		//we still want the key to be found in case it is alone on the line
		vector<string>& contents = strings[key];
		
		while (!source.eof()) {
			//On 'IX machines, sometimes that extra windows return thing will 
			//prevent us from seeing the end of line on time....gotta clean out
			//the old value
			word = "";					
			source>>word;
			if (word.length() > 0) {
				contents.push_back(word);
			}
		}
	}
	return valid;
}
		
inline
FileToArray::FileToArray() {}

inline
FileToArray::~FileToArray() {}

inline
bool FileToArray::ParseLine(const char *line, uint val) {
	if (strlen(line) > 0) {
		strings.push_back(string(line));
		return true;
	}
	return false;
}

inline
StringArray FileToArray::OpenFile(const char *filename, char cmtChar) {
	LineParser lp(cmtChar);
	lp.Parse(filename, this, false);
	return strings;
}

inline
bool BasicLineCounter::ParseLine(const char *line, uint val) {
	if (strlen(line) > 0 && line[0] != '/' && line[0] != '#')
		return true;
	else
		return false;
}

inline
uint FileToMap::Parse(const char *filename) {
	LineParser lp;
	keys.clear();
	uint linesKept = lp.Parse(filename, this);
//	if (keepComments)
//		comments = lp.comments;
	return linesKept;
}

inline
uint FileToMapExplicit::Parse(const char *filename) {
	LineParser lp;
	InitKey("PROJECT",				filename,							"#The name of the project");

	uint linesKept = lp.Parse(filename, this);
//	if (keepComments)
//		comments = lp.comments;
	return linesKept;
}
inline
void FileToMapExplicit::InitKey(const char *key, const char *defValue, const char *comment) {
	string line = string(key) + string(" ") + string(defValue);

	//Add the key, if it's not already there
	if (strings.find(key) == strings.end())
		keys.push_back(key);

	ParseLine(line.c_str());
	defaultValue[key]=true;
	varComments[key] = comment;
}

inline
void FileToMap::Write(ostream& os) {
	vector<string>::iterator kItr = keys.begin();
	vector<string>::iterator kEnd = keys.end();
	map<string, string>::iterator commentEnd = varComments.end();
	while (kItr != kEnd) {
		if (varComments.find(*kItr) != commentEnd) 
			os<<"\n"<<comment<<" "<<varComments[*kItr]<<"\n";
		vector<string> strings;
		GetLines((*kItr).c_str(), strings);
		vector<string>::iterator lItr = strings.begin();
		vector<string>::iterator lEnd = strings.end();
		if (lItr==lEnd) 
			os<<"#"<<*kItr<<"\n";
		while (lItr!=lEnd) {
			os<<*kItr<<"\t"<<*lItr<<"\n";
			lItr++;
		}
		kItr++;
	}
}

inline
void FileToMapExplicit::AppendValue(const char *key, const char *val) {
	vector<string>& list = strings[key];
	if (defaultValue[key])
		list.clear();
	defaultValue[key] = false;
	list.push_back(val);
}

inline
bool FileToMapExplicit::ParseLine(const char *line, uint val) {
	stringstream source(line);

	bool valid = strlen(line) > 0 && line[0] != '\r' && line[0] != '\n';
	if (line[0] == comment) {
		comments+=string(line) + string("\n");
		return false;
	}
	if (valid) {
		string key, word;

		source>>key;
	
		if (find(keys.begin(), keys.end(), key) == keys.end())  {
			cout<<"!! Configuration Parser Error. Unknown Key, '"<<key<<"'.\n>>"<<line<<"\n";
			return false;
		}

		if (comments.length() > 0)
			varComments[key]=varComments[key]+comments;
		comments="";
		
		//We don't really care if contents is empty, 
		//we still want the key to be found in case it is alone on the line
		vector<string>& contents = strings[key];

		if (defaultValue[key])
			contents.clear();
		defaultValue[key]=false;

		while (!source.eof()) {
			//On 'IX machines, sometimes that extra windows return thing will 
			//prevent us from seeing the end of line on time....gotta clean out
			//the old value
			word = "";					
			source>>word;
			if (word.length() > 0) {
				contents.push_back(word);
			}
		}


	}
	return valid;
}
	
}

#endif
