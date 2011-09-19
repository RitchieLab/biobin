/* 
 * File:   taskfilegeneration.h
 * Author: torstees
 *
 * Created on June 30, 2011, 9:50 AM
 */

#ifndef TASKFILEGENERATION_H
#define	TASKFILEGENERATION_H

#include "biofilter/task.h"
#include "binapplication.h"

namespace BioBin {
namespace Task {


class GenerateFiles : public Biofilter::Task::Task {
public:
	GenerateFiles();
	virtual ~GenerateFiles();

	virtual void Init(Biofilter::Application* app);

	void ExecuteTask();

	static bool WriteBinData;
	static bool WriteGenotypeData;
	static std::string OutputDelimeter;
protected:
	BinApplication* app;
};


inline
GenerateFiles::GenerateFiles() : Task(3), app(NULL) {}

inline
GenerateFiles::~GenerateFiles() { }

inline
void GenerateFiles::Init(Biofilter::Application* app) {
	this->app = (BinApplication*)app;
}


}
}

#endif	/* TASKFILEGENERATION_H */

