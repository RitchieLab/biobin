//
// C++ Interface: exception
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYEXCEPTION_H
#define UTILITYEXCEPTION_H

#include <sys/types.h>
#include <string>
#include <iostream>
#include "strings.h"

namespace Utility {

namespace Exception {

using namespace std;
/**
Base class for all exceptions

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Base{
public:
    Base() { //cerr<<"Exception created\n"; 
	}

    virtual ~Base() { }

	virtual string GetErrorMessage() const { return "";}
	friend ostream & operator << (ostream& os, const Base& e) { os<<e.GetErrorMessage(); return os; }
private:

};

class General : public Base {
public:
	General(const char *message) : message(message) { }
	General(string& message) : message(message) { }
	General() { }
	string GetErrorMessage() const { return message; }
protected:
	string message;
};

class IndexOutOfRange : public General {
public:
	IndexOutOfRange(const char *variableName, uint size = 0, uint value = 0) : General("") {
		message = "Attempt to access " + std::string(variableName) + " beyond it's bounds. ";
		if (size > 0 || value > 0)
			message += " Size: " + ToString<uint>(size) + "\tRequested: " + ToString<uint>(value);
		message += "\n";
	}
};

class MissingKey : public General {
public:
	MissingKey(const char *variableName, const char *key, uint size) {
		message = "Failed attempt to lookup " + std::string(key) + " from " + std::string(variableName);
		if (size > 0)
			message += " which contains " + ToString(size) + " entries";
		message += ".\n";
	}
};

class FileIO : public General {
public:
	FileIO(const char *filename, const char *message) : General(message), filename(filename) { }
	FileIO(const char *filename) : filename(filename) { }
	string filename;
};

class FileNotFound: public FileIO {
public:
	FileNotFound(const char *filename) : FileIO(filename) { 
		message = "The file, '" + this->filename + "', couldn't be found. ";
	}
};

class FileNotWritable: public FileIO {
public:
	FileNotWritable(const char *filename) : FileIO(filename) {
		message = "The file, " + this->filename + ", couldn't be written to. ";
	}
};

class BadFileVersion: public FileIO {
public:
	BadFileVersion(const char *filename)  : FileIO(filename) {
		message = "The file, " + this->filename + ", is incompatible with the current application. Please verify that you have selected the correct file.";
	}
};

}
}
#endif
