#define BOOST_TEST_MODULE order_of_operations
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

EXPRTEST(pemdas1, "2*3+4*5",   2*3+4*5)
EXPRTEST(pemdas2, "2*(3+4)*5", 2*(3+4)*5)
EXPRTEST(pemdas3, "2**3+4",    std::pow(2,3)+4)
EXPRTEST(pemdas4, "2**2**-3",  std::pow(2.,std::pow(2.,-3.)))

EXPRTEST(rel1 , "1 == 1",   1 == 1)
EXPRTEST(rel2 , "0 == 1",   0 == 1)
EXPRTEST(rel3 , "1 == 0",   1 == 0)

EXPRTEST(rel4 , "1 != 1",   1 != 1)
EXPRTEST(rel5 , "0 != 1",   0 != 1)
EXPRTEST(rel6 , "1 != 0",   1 != 0)

EXPRTEST(rel7 , "1 >= 1",   1 >= 1)
EXPRTEST(rel8 , "0 >= 1",   0 >= 1)
EXPRTEST(rel9 , "1 >= 0",   1 >= 0)

EXPRTEST(rel10, "1 <= 1",   1 <= 1)
EXPRTEST(rel11, "0 <= 1",   0 <= 1)
EXPRTEST(rel12, "1 <= 0",   1 <= 0)

EXPRTEST(rel13, "1 > 1" ,   1 > 1 )
EXPRTEST(rel14, "0 > 1" ,   0 > 1 )
EXPRTEST(rel15, "1 > 0" ,   1 > 0 )

EXPRTEST(rel16, "1 < 1" ,   1 < 1 )
EXPRTEST(rel17, "0 < 1" ,   0 < 1 )
EXPRTEST(rel18, "1 < 0" ,   1 < 0 )

EXPRTEST(rel19, "1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1", 1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1)
EXPRTEST(rel20, "0 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1", 0 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1)
EXPRTEST(rel21, "1 && 0 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1", 1 && 0 == 1 != 1 >= 1 <= 1 > 1 < 1 || 1)
EXPRTEST(rel22, "1 && 1 == 0 != 1 >= 1 <= 1 > 1 < 1 || 1", 1 && 1 == 0 != 1 >= 1 <= 1 > 1 < 1 || 1)
EXPRTEST(rel23, "1 && 1 == 1 != 0 >= 1 <= 1 > 1 < 1 || 1", 1 && 1 == 1 != 0 >= 1 <= 1 > 1 < 1 || 1)
EXPRTEST(rel24, "1 && 1 == 1 != 1 >= 0 <= 1 > 1 < 1 || 1", 1 && 1 == 1 != 1 >= 0 <= 1 > 1 < 1 || 1)
EXPRTEST(rel25, "1 && 1 == 1 != 1 >= 1 <= 0 > 1 < 1 || 1", 1 && 1 == 1 != 1 >= 1 <= 0 > 1 < 1 || 1)
EXPRTEST(rel26, "1 && 1 == 1 != 1 >= 1 <= 1 > 0 < 1 || 1", 1 && 1 == 1 != 1 >= 1 <= 1 > 0 < 1 || 1)
EXPRTEST(rel27, "1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 0 || 1", 1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 0 || 1)
EXPRTEST(rel28, "1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 0", 1 && 1 == 1 != 1 >= 1 <= 1 > 1 < 1 || 0)
