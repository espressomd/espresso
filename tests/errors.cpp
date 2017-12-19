#define BOOST_TEST_MODULE errors
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "matheval.hpp"

BOOST_AUTO_TEST_CASE(parse_failure)
{
    std::string const s = "#";
    BOOST_CHECK_THROW(matheval::parse<double>(s, {}),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(expectation_failure)
{
    std::string const s = "(";
    BOOST_CHECK_THROW(matheval::parse<double>(s, {}),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(unknown_variable)
{
    std::string const s = "x";
    BOOST_CHECK_THROW(matheval::parse<double>(s, {}),
                      std::invalid_argument);
}
