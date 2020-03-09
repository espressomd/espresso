#define BOOST_TEST_MODULE algorithms_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <limits>

#include <observables/algorithms.hpp>

using namespace Observables::Algorithms;

BOOST_AUTO_TEST_CASE(algorithms) {
  auto identity = [](auto const &p) { return p; };
  std::vector<int> values{1, 2, 3, 4};
  {
    auto one = [](auto const &p) { return 1; };
    auto const res = WeightedAverage()(values, identity, one);
    BOOST_CHECK(res == std::accumulate(values.begin(), values.end(), 0) /
                           values.size());
  }
  {
    auto plus_one = [](auto const &p) { return p + 1; };
    auto const res = WeightedAverage()(values, identity, plus_one);
    BOOST_CHECK(res == (1 * 2 + 2 * 3 + 3 * 4 + 4 * 5) / 14);
    auto const res2 = WeightedSum()(values, identity, plus_one);
    BOOST_CHECK(res2 == (1 * 2 + 2 * 3 + 3 * 4 + 4 * 5));
  }
  {
    auto const res = Average()(values, identity);
    BOOST_CHECK(res == std::accumulate(values.begin(), values.end(), 0) /
                           values.size());
  }
}
