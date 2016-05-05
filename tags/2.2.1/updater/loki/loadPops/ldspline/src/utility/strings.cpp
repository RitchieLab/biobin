#include "strings.h"
#include <sstream>
#include <iomanip>
#include <iostream>

#include <sys/stat.h>
#include <sstream>
#include "exception.h"
#include <fstream>
#include <cctype>
#include <algorithm>
#include <cstring>

namespace LDUtility {

#ifdef WIN32 
const char *seperators ="\\/";
#else
const char *seperators = "\\/";
#endif


struct comma_facet : public std::numpunct<char>
{ protected: string do_grouping() const { return "\003" ; } };



string ThousandsFormat(double d) {
	stringstream ss;
	locale l = ss.getloc();
	ss.imbue( locale(l, new comma_facet));
	ss<<setiosflags(ios::fixed)<<d;
	return ss.str();
}

string StripExtension(const char *originalFilename) {
	string filename(originalFilename);

	size_t firstDot = filename.find_first_of(".");
	if (firstDot != string::npos)
		filename.erase(firstDot, string::npos);
	return filename;
}
/**
 * @brief Extracts the filename (right of the rightmost "/" and left of the leftmost "."
 */
string ExtractBaseFilename(const char *originalFilename) {
	string filename(originalFilename);

	size_t lastSlash = filename.find_last_of(seperators);
	if (lastSlash != string::npos)
		filename.erase(0, lastSlash + 1);

	size_t lastDot = filename.find_last_of(".");
	if (lastDot != string::npos)
		filename.erase(lastDot);

	return filename;
}

uint ExtractRsNumber(const char *rsid) {
	string rs(rsid);
	if (rs.find("r") != string::npos || rs.find("R") != string::npos)
			rs.erase(0,2);
	return atoi(rs.c_str());
}

string ExtractFilename(const char *maybeHasFullPath) {
	string filename(maybeHasFullPath);
	
	size_t lastSlash = filename.find_last_of(seperators);
	if (lastSlash != string::npos)
		filename.erase(0, lastSlash + 1);

	return filename;
}

string EscapeSpaces(const char *filename ){ 
	stringstream ss(filename);
	
	stringstream newName;

	string word;
	ss>>word;
	newName<<"\"";
/*	newName<<word;
	while (!ss.eof() ) {
		ss>>word;
		newName<<"\\ "<<word;
	}
*/	newName<<filename<<"\"";
		
	return newName.str();
}

string ExtractExtension(const char *filename) {
	string extension(filename);
	size_t lastDot = extension.find_last_of(".");
	if (lastDot != string::npos)
		extension.erase(0, lastDot+1);
	return extension;
}

string ParseFilename(istream &s, const char *desc) {
	string filename;
	s>>filename; 

	size_t lastPos = filename.length() - 1;
	if (filename[0] == '"' && filename[lastPos] !='"') {
		//Parsing filename with spaces
		bool doContinue = true;
		while (doContinue && !s.eof()) {
			string word;
			s>>word;
			filename+=" " + word;
			doContinue = word[word.length() -1]!='"';
		}
		if (doContinue) {
			stringstream ss;
			ss<<"Unterminated quote was encountered in the configuration "<<desc<<"\n"<<filename<<"\n";
			throw Exception::General(ss.str().c_str());
		}
	}
	return StripQuotes(filename.c_str());
}

/**
 * @brief attempts to break the filename into 3 components
 */
void SplitIntoComponents(const char *originalFilename, string &filepath, string &basefilename, string &ext) {
	basefilename = ext = filepath = originalFilename;

	size_t lastDot = basefilename.find_last_of(".");
	if (lastDot != string::npos) {
		basefilename.erase(lastDot);
		ext.erase(0, lastDot + 1);
	}
	else
		ext = "";

	size_t lastSlash = basefilename.find_last_of(seperators);
	if (lastSlash != string::npos) {
		filepath.erase(lastSlash+1);
		basefilename.erase(0, lastSlash + 1);
	}
	else {
		filepath="";
		lastSlash = 0;
	}
}

// Quick way to remove whitespace from strings (windows /r for the most part from end of lines)
// http://www.cplusplus.com/reference/string/string/find_last_not_of/
string StripTrailingWhitespace(const char *word){ 
	string whitespace("\t\f\v\n\r ");
	string str = word;
	size_t found=str.find_last_not_of(whitespace);
	if (found!=string::npos)
		return str.erase(found+1);
	else 
		return "";
}
	
string StripQuotes(const char *Filename) {
	string filename(Filename);
	if (filename[0] == '"')
		filename.erase(0, 1);
	if (filename[filename.length() - 1] == '"')
		filename.erase(filename.length() - 1);
	return filename;
}

/**
 * @brief Extracts tokens seperated according the characters passed as sep
 */
uint TokenizeString(const char *origString, vector<string>& tokens, const char *sep) {
	string filename(origString);
	size_t strLength = filename.length();
	char *tok = new char[strLength + 1];
	size_t firstToken = 0;
	size_t nextToken = filename.find_first_of(sep);

	tokens.clear();

	//std::cout<<"Parsing string: "<<origString<<" by ("<<sep<<")\n";

	while (firstToken != string::npos) {
		int tokLength = nextToken - firstToken;
//		if (nextToken < ) {
			if (nextToken != string::npos || nextToken < strLength) {
				strncpy(tok, origString + firstToken, tokLength);
				nextToken++;
			}
			else {
				strcpy(tok, origString + firstToken);
				tokLength = strLength - firstToken;
			}
	
			tok[tokLength] = '\0';
			if (tokLength > 0) {
//				std::cout<<"\t"<<count++<<" "<<tok<<"\n";
				tokens.push_back(tok);
			}
//		}
		firstToken = nextToken;
		nextToken = filename.find_first_of(sep, firstToken);

	}
	return tokens.size();
}


uint CountColumns(const char *line) {
	stringstream ss(line);
	string vals;

	uint count = 0;

	while (!ss.eof()) {
		vals = "xvx";
		ss>>vals;	
		
		if (vals != "xvx")
			count++;
	}
	
	return count;
}

bool FileExists(const char *filename) {

	ifstream testFile(filename);
	return testFile.is_open();

	struct stat stFileInfo;
/*	cout<<"Checking for the existance of "<<filename;
	if (stat(filename, &stFileInfo) == 0)
		cout<<"Yep, it's there!\n";
	else 
		cout<<"Nope. Can't find anything of that sort!\n";
 */
	string fixedFilename = EscapeSpaces(filename);
	if (stat(fixedFilename.c_str(), &stFileInfo))
		cout<<"Can't find file, "<<fixedFilename<<"\n";
	return stat(fixedFilename.c_str(),&stFileInfo) == 0;
}


string FileToString(const char *filename, const char *sep) {
	string contents;
	ifstream file(filename);
	while (file.good() && !file.eof()) {
		string word = "";
		file>>word;
		if (word.length() > 0) {
			if (contents.length() > 0)
				contents += sep;
			contents+=word;
		}
	}
	return contents;
}

vector<string> FileToVector(const char *filename) {
	vector<string> contents;
	ifstream file(filename);
	while (file.good() && !file.eof()) {
		string word = "";
		file>>word;
		if (word.length() > 0) 
			contents.push_back(word);
	}
	return contents;
}
string ToString(int value) {
	char s[64];
	sprintf(s, "%d", value);
	return s;
}

string ToString(uint value) {
	char s[64];
	sprintf(s, "%d", value);
	return s;
}

string ToString(double value, int prec) {
	stringstream s;
	
	if (prec>0)
		s<<setprecision(prec);
	s<<value;
	return s.str();
}

string ToUpper(const char *mixedcase) {
	string allupper = mixedcase; 
	transform(allupper.begin(), allupper.end(), allupper.begin(), (int(*)(int))toupper);
	return allupper;
}

string IntToChrom(int chrom) {
	static string labels[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};
	return labels[chrom-1];
}

int ChromToInt(const char *chrom) {
	if (strcmp("X", chrom)==0 || strcmp("X ", chrom)==0 || strcmp(" X", chrom)==0)
		return 23;
	else if (strcmp("Y", chrom)==0 || strcmp("Y ", chrom)==0 || strcmp(" Y", chrom)==0)
		return 24;
	else if (strcmp("MT", chrom)==0)
		return 25;
	else
		return atoi(chrom);
}

}

