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

#include <string>
#include <iostream>

namespace LDUtility {

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

class FileIO : public General {
public:
	FileIO(const char *filename, const char *message) : General(message), filename(filename) { }
	FileIO(const char *filename) : filename(filename) { }
	string filename;
};

class FileNotFound: public FileIO {
public:
	FileNotFound(const char *filename) : FileIO(filename) { 
		message = "The file, " + this->filename + ", couldn't be found. ";
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
