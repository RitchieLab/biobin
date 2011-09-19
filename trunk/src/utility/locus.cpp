/* 
 * File:   locus.cpp
 * Author: torstees
 * 
 * Created on July 1, 2011, 3:45 PM
 */

#include "locus.h"

namespace Utility {
//uint Locus::MaxChromosomes = 26;

// These should be set up based on the presence of certain words in the various 
// roles from the database.
Utility::IdCollection Locus::ExonRoles;
Utility::IdCollection Locus::IntronRoles;
Utility::IdCollection Locus::RegulatoryRoles;

#ifdef TEST_APP

// Yep, not much here...not enough time to write tests 
TEST(RSTest, ChromConversion) {
	EXPECT_EQ(1, ChromToInt("1 "));
	EXPECT_EQ(1, ChromToInt("Chr1"));
	EXPECT_EQ(1, ChromToInt("chr1 "));
	EXPECT_EQ(2, ChromToInt("CHR2 "));
	EXPECT_EQ(26, ChromToInt("CHRMT"));
	EXPECT_EQ(22, ChromToInt("22"));
	EXPECT_EQ(26, ChromToInt("MT"));
	EXPECT_EQ(24, ChromToInt("Y "));
	EXPECT_EQ(23, ChromToInt(" X"));
}

#endif 

}