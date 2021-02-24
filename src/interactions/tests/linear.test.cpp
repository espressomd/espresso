#define BOOST_TEST_MODULE interactions linear
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interactions/linear.hpp"

constexpr auto tol = std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_CASE(linear) {
  Interactions::Linear<double> lin{3.5, 1.0};
  BOOST_CHECK_SMALL(lin.energy(0.0) - 1, tol);
  BOOST_CHECK_SMALL(lin.energy(1.0) - (-3.5 + 1.0), tol);
}