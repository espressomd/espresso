#define BOOST_TEST_MODULE order_of_operations
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

EXPRTEST(pemdas1, "2*3+4*5",   2*3+4*5)
EXPRTEST(pemdas2, "2*(3+4)*5", 2*(3+4)*5)
EXPRTEST(pemdas3, "2**3+4",    std::pow(2,3)+4)
EXPRTEST(pemdas4, "2**2**-3",  std::pow(2.,std::pow(2.,-3.)))
