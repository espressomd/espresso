#define BOOST_TEST_MODULE acumulator test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Accumulator.hpp"

BOOST_AUTO_TEST_CASE(accumulator) {
  auto acc = Utils::Accumulator(4);
  auto test_data1 = std::vector<double>{{0.0, 1.0, 2.0, 3.0}};
  auto test_data2 = std::vector<double>{{1.5, 3.5, 3.5, 4.5}};
  acc(test_data1);
  BOOST_CHECK(acc.get_mean() == test_data1);
  BOOST_CHECK(acc.get_variance() ==
              std::vector<double>(4, std::numeric_limits<double>::max()));
  acc(test_data2);
  BOOST_CHECK(
      (acc.get_mean() == std::vector<double>{{0.75, 2.25, 2.75, 3.75}}));
  BOOST_CHECK((acc.get_variance() ==
               std::vector<double>{{1.125, 3.125, 1.125, 1.125}}));
  BOOST_CHECK(
      (acc.get_std_error() == std::vector<double>{{0.75, 1.25, 0.75, 0.75}}));
}
