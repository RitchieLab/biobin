/* 
 * File:   taskfilegeneration.h
 * Author: torstees
 *
 * Created on June 30, 2011, 9:50 AM
 */

#ifndef TASKFILEGENERATION_H
#define	TASKFILEGENERATION_H

#include "task.h"
#include "binapplication.h"

namespace BioBin {

class Application;

namespace Task {


class GenerateFiles : public Task {
public:
	GenerateFiles();
	virtual ~GenerateFiles();

	virtual void Init(Application* app);

	void ExecuteTask();

	static bool WriteBinData;
	static bool WriteGenotypeData;
	static std::string OutputDelimeter;
protected:
	BinApplication* app;
};





}
}

#endif	/* TASKFILEGENERATION_H */

