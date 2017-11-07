#include <string>
#include "matheval.hpp"

#define EXPRTEST(casename, expr, expected)                             \
BOOST_AUTO_TEST_CASE( casename )                                       \
{                                                                      \
    std::string const s = expr;                                        \
    std::map<std::string, double> st;                                  \
    double result;                                                     \
    BOOST_CHECK_NO_THROW(result = matheval::parse<double>(s, st));     \
    BOOST_CHECK_CLOSE(result, (expected),                              \
                      std::numeric_limits<double>::epsilon());         \
}

#define SYMEXPRTEST(casename, expr, st, expected)                      \
BOOST_AUTO_TEST_CASE( casename )                                       \
{                                                                      \
    std::string const s = expr;                                        \
    double result;                                                     \
    BOOST_CHECK_NO_THROW(result = matheval::parse<double>(s, st));     \
    BOOST_CHECK_CLOSE(result, (expected),                              \
                      std::numeric_limits<double>::epsilon());         \
}
