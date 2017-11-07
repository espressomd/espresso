#define BOOST_TEST_MODULE literals
#include <boost/test/included/unit_test.hpp>
#include "exprtest.hpp"

EXPRTEST(literal1, "1.234",    1.234)
EXPRTEST(literal2, "4.2e2",    420)
EXPRTEST(literal3, "5e-01",    0.5)
EXPRTEST(literal4, "-3",      -3)
EXPRTEST(literal5, "pi",       boost::math::constants::pi<double>())
EXPRTEST(literal6, "epsilon",  std::numeric_limits<double>::epsilon())
EXPRTEST(literal7, "digits",   std::numeric_limits<double>::digits)
EXPRTEST(literal8, "digits10", std::numeric_limits<double>::digits10)
EXPRTEST(literal9, "e",        boost::math::constants::e<double>())
