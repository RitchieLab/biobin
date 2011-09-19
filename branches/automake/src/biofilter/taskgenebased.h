/* 
 * File:   taskgenebased.h
 * Author: torstees
 *
 * Created on April 15, 2011, 11:07 AM
 */

#ifndef TASKGENEBASED_H
#define	TASKGENEBASED_H

#include "utility/types.h"
#include "task.h"
#include <map>

namespace Biofilter {

namespace Task {

class GeneBased : public Task {
public:
	GeneBased();


	virtual ~GeneBased() {}

	virtual void Init(Application* app);

protected:
	std::multimap<uint, uint> *geneLookup;
};

inline
GeneBased::GeneBased() : Task(2),
		geneLookup(NULL) {}


inline
void GeneBased::Init(Application* app) {
	regions			= app->GetRegions();
	snps				= app->GetDataset();
	geneLookup		= app->GetGeneLookup();
	filename			= app->AddReport(GetFileSuffix().c_str(), GetFileExtension().c_str(), GetFileDescription().c_str());
}

}
}
#endif	/* TASKGENEBASED_H */

