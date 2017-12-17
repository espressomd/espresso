#define BOOST_TEST_MODULE literals
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "exprtest.hpp"

using namespace boost::math::constants;

EXPRTEST(literal1, "1.234"  ,  1.234)
EXPRTEST(literal2, "4.2e2"  ,  420)
EXPRTEST(literal3, "5e-01"  ,  0.5)
EXPRTEST(literal4, "-3"     , -3)
EXPRTEST(literal5, "pi"     ,  pi<double>())
EXPRTEST(literal6, "epsilon",  std::numeric_limits<double>::epsilon())
EXPRTEST(literal7, "phi"    ,  phi<double>())
EXPRTEST(literal8, "e"      ,  e<double>())
