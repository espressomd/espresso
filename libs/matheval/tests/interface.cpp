#define BOOST_TEST_MODULE interface
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

#include "matheval.hpp"

BOOST_AUTO_TEST_CASE(integration1) {
    std::string expr = "pow(x/2 + sqrt(x**2/4 + y**3/24), 1/3)";
    double x = 2.0, y = -1.0;

    matheval::Parser parser;
    BOOST_CHECK_NO_THROW(parser.parse(expr));

    std::map<std::string, double> symbol_table;
    symbol_table.insert(std::make_pair("x", x));
    symbol_table.insert(std::make_pair("y", y));

    double result = 0;
    BOOST_CHECK_NO_THROW(result = parser.evaluate(symbol_table));

    double expected = std::pow(
        x / 2. + std::sqrt(std::pow(x, 2.) / 4. + std::pow(y, 3.) / 24.),
        1. / 3.);
    BOOST_CHECK_CLOSE_FRACTION(result, expected,
                               std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(integration2) {
    std::string expr = "(";

    matheval::Parser parser;

    // Parsing should fail
    BOOST_CHECK_THROW(parser.parse(expr), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(integration3) {
    std::string expr = "x";

    matheval::Parser parser;

    BOOST_CHECK_NO_THROW(parser.parse(expr));

    // Evaluating should fail
    BOOST_CHECK_THROW(parser.evaluate(),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(integration4) {
    std::string expr = "1 + 1";

    matheval::Parser parser;

    BOOST_CHECK_NO_THROW(parser.parse(expr));

    BOOST_CHECK_NO_THROW(parser.optimize());
}
