#define BOOST_TEST_MODULE arithmetic
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

EXPRTEST(basicop1, " 2 +\t3\n",  5)       // Whitespace ignored
EXPRTEST(basicop2, " 2 -\t3\n", -1)
EXPRTEST(basicop3, " 2 *\t3\n",  6)
EXPRTEST(basicop4, " 2 /\t3\n",  2./3.)   // Double division
EXPRTEST(basicop5, " 2 ** 3\n",  8)
