#include "types.h"
#include <stdlib.h>
#ifdef WIN32

float drand48 ()
{
	return (float)(((double) (rand ( ))) / ((double) RAND_MAX)) ;
}
/*
float nanf(const char *tagp) {
	if (strcmp(tagp, "(charsequence)") == 0)
		return strtod("NAN(charsequence)", NULL);

	//I'm not sure what else we might return!
	return strtod("NAN(charsequence)", NULL);
}
*/	
#endif //WIN32
