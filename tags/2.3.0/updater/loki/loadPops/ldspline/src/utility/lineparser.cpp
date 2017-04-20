//
// C++ Implementation: lineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineparser.h"
#include "exception.h"
namespace LDUtility {

uint LineParser::Parse(const char *filename, AsciiParser *parser, bool keepComments) {
	char line[MAX_LINE_LENGTH];
	uint linesParsed = 0;							//Keep up with what we've seen
	uint lineCount	= 0;							//Used to know how many were skipped

	ifstream file(filename, ios_base::in);			//open file for reading
	if (file.is_open()) {
		while (!file.eof()) {

			file.getline(line, MAX_LINE_LENGTH);
			//Let's ignore anything after a the comment character
	 		if (cmtChar) {
				char *cmtSection = strchr(line, cmtChar);
				if (cmtSection) {
					if (keepComments)
						comments+=string(line)+string("\n");
					(*cmtSection) = '\n';
					memset((void*)(cmtSection + 1), ' ', strlen(line)-(cmtSection-line) - 1);
				}
			}
			if (parser->ParseLine(line, linesParsed))
				lineCount++;
			linesParsed++;
		}
	}
	else {
		cout<<"File `"<<filename<<"` wasn't opened\n";
		throw Exception::FileNotFound(filename);
	}
		
	return lineCount;
}



}//Utility
