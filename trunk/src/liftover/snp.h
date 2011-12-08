/* 
 * File:   snp.h
 * Author: torstees
 *
 * Created on May 6, 2011, 1:28 PM
 */

#ifndef SNP_H
#define	SNP_H

#include <sys/types.h>
#include <vector>
#include "utility/strings.h"
#include "utility/locus.h"

namespace LiftOver {

typedef Utility::Locus SNP;
typedef std::vector<Utility::Locus> SnpArray;


}
#endif	/* SNP_H */

