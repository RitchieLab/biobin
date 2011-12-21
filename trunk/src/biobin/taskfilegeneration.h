/* 
 * File:   taskfilegeneration.h
 * Author: torstees
 *
 * Created on June 30, 2011, 9:50 AM
 */

#ifndef TASKFILEGENERATION_H
#define	TASKFILEGENERATION_H

#include "task.h"

namespace BioBin {

class BinApplication;

namespace Task {


class GenerateFiles : public Task {
public:
	GenerateFiles(BinApplication* app);
	virtual ~GenerateFiles(){}

	virtual void ExecuteTask();

	static bool WriteBinData;
	static bool WriteGenotypeData;
	static bool WriteLociData;
	static std::string OutputDelimiter;
protected:
	BinApplication* _app;
};





}
}

#endif	/* TASKFILEGENERATION_H */

