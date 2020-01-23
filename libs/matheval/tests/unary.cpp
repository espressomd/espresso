#define BOOST_TEST_MODULE unary
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

EXPRTEST(unary1, "-(2)",  -2)
EXPRTEST(unary2, "-(-2)",  2)
EXPRTEST(unary3, "+(-2)", -2)
EXPRTEST(unary4, "+(+2)",  2)
EXPRTEST(unary5, "!(1)",   0)
EXPRTEST(unary6, "!(0)",   1)
