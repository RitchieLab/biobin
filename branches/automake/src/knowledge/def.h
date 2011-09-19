/* 
 * File:   def.h
 * Author: torstees
 *
 * Common Definitions to make it easier to find them
 * Created on March 7, 2011, 4:57 PM
 */

#ifndef DEF_H
#define	DEF_H

#include "utility/types.h"

namespace Knowledge {

#define MAX_RS_LENGTH 12		///< Used pad the text file
#define MAX_II_LENGTH 7			///< Length of II field (snpsnp text)

extern bool BinaryArchive;		///< Used to control output type for the various files

namespace ModelGenerationMode {
	typedef uint Type;
	const Type ALL_MODELS				= 1;
	const Type GROUP_LEVEL				= 2;
	const Type DD_ONLY					= 3;
	const Type UNRECOGNIZED				= (uint)-1;

	Type ConvertType(const char* name);

};

namespace MetaGroup {
	typedef uint Type;
	const uint DiseaseIndependent		= 1;
	const uint DiseaseDependent		= 2;
	const uint SnpCollection			= 3;
	const uint GeneCollection			= 4;
	const uint MetaGroupCount			= 5;
	const uint UNRECOGNIZED				= (uint)-1;

	Type ConvertType(const char* name);
	std::string ConvertType(MetaGroup::Type t);
};

inline
ModelGenerationMode::Type ModelGenerationMode::ConvertType(const char* name) {
	if (strcmp(name, "ALL_MODELS")==0)
		return ALL_MODELS;
	if (strcmp(name, "GROUP_LEVEL")==0)
		return GROUP_LEVEL;
	if (strcmp(name, "DD_ONLY")==0)
		return DD_ONLY;
	return UNRECOGNIZED;
}

inline
std::string MetaGroup::ConvertType(MetaGroup::Type t) {
	switch(t) {
	case DiseaseIndependent:
		return "DISEASE_INDEPENDENT";
	case DiseaseDependent:
		return "DISEASE_DEPENDENT";
	case SnpCollection:
		return "SNP_COLLECTION";
	case GeneCollection:
		return "GENE_COLLECTION";
	case UNRECOGNIZED:
		return "Unrecognized Type";
	}
	return "Unrecognized Type";
}

inline
MetaGroup::Type MetaGroup::ConvertType(const char* name) {
		if (strcmp(name, "DISEASE_INDEPENDENT")==0)
			return DiseaseIndependent;
		if (strcmp(name, "DISEASE_DEPENDENT")==0)
			return DiseaseDependent;
		if (strcmp(name, "SNP_COLLECTION")==0)
			return SnpCollection;
		if (strcmp(name, "GENE_COLLECTION")==0)
			return GeneCollection;
		return UNRECOGNIZED;
	}
}
#endif	/* DEF_H */

