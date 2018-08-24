#include <algorithm>
#include <vector>
#include <iterator>
#include <numeric>

#define BOOST_TEST_MODULE SkipIterator test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/utils/SkipIterator.hpp"

BOOST_AUTO_TEST_CASE(skip) {
  using Utils::make_skip_iterator;
  
  std::vector<int> data(50);
  std::iota(data.begin(), data.end(), 1);

  auto multiple_of_3 = [](int i) -> bool { return (i % 3) == 0; };

  auto begin = make_skip_iterator(data.begin(), data.end(), multiple_of_3);
  auto end = make_skip_iterator(data.end(), data.end(), multiple_of_3);

  for(; begin != end; ++begin) {
    BOOST_CHECK(!multiple_of_3(*begin));
  }
}

BOOST_AUTO_TEST_CASE(empty) {
  using Utils::make_skip_iterator;
  
  std::vector<int> data(50);
  std::iota(data.begin(), data.end(), 1);

  auto always_true = [](int) -> bool { return true; };

  auto begin = make_skip_iterator(data.begin(), data.end(), always_true);
  auto end = make_skip_iterator(data.end(), data.end(), always_true);

  BOOST_CHECK(begin == end);
}

BOOST_AUTO_TEST_CASE(gap) {
  using Utils::make_skip_iterator;
  
  std::vector<int> values;

  for(int i = 1; i <= 10; i++) {
    values.push_back(i);
    for(int j = 0; j < i; j++) {
      values.push_back(0);
    }
  }

  auto zero = [](int i) -> bool { return i == 0; };

  auto begin = make_skip_iterator(values.begin(), values.end(), zero);
  auto end = make_skip_iterator(values.end(), values.end(), zero);

  BOOST_CHECK(std::distance(begin, end) == 10);
}
