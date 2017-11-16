/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define BOOST_TEST_MODULE Utils::Cache test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/make_unique.hpp"
using Utils::make_unique;

#include "utils/Cache.hpp"
using Utils::Cache;
using Utils::make_cache;

BOOST_AUTO_TEST_CASE(types) {
  auto null = [](int) { return make_unique<const char>(0); };
  using cache_type = Cache<int, char, decltype(null)>;
  static_assert(std::is_same<cache_type::key_type, int>::value, "");
  static_assert(std::is_same<cache_type::value_type, const char *>::value, "");
}

BOOST_AUTO_TEST_CASE(get_value) {
  {
    auto next = [](int i) { return make_unique<const char>(i + 1); };

    Cache<int, char, decltype(next)> cache{next};
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(42) == 43);
  }
  {
    auto null_getter = [](int) { return std::unique_ptr<const char>{}; };
    Cache<int, char, decltype(null_getter)> cache{null_getter};
    BOOST_CHECK(!cache.get(123));
    BOOST_CHECK(!cache.get(11));
  }
}

BOOST_AUTO_TEST_CASE(has) {
  {
    auto next = [](int i) { return make_unique<const char>(i + 1); };

    Cache<int, char, decltype(next)> cache{next};
    cache.get(41);
    cache.get(1);

    BOOST_CHECK(cache.has(41));
    BOOST_CHECK(cache.has(1));
    BOOST_CHECK(!cache.has(2));
    BOOST_CHECK(cache.has(1));
  }
}

BOOST_AUTO_TEST_CASE(invalidate) {
  {
    auto next = [](int i) { return make_unique<const char>(i + 1); };

    Cache<int, char, decltype(next)> cache{next};
    cache.get(41);
    cache.get(1);

    cache.invalidate();

    BOOST_CHECK(!cache.has(41));
    BOOST_CHECK(!cache.has(1));
    BOOST_CHECK(!cache.has(2));
  }
}

BOOST_AUTO_TEST_CASE(caching) {
  {
    auto count = 0;
    auto counter = [&count](int i) {
      count++;
      return make_unique<const char>(i + 1);
    };

    Cache<int, char, decltype(counter)> cache{counter};
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);

    BOOST_CHECK(1 == count);
    cache.invalidate();

    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);
    BOOST_CHECK(*cache.get(41) == 42);

    BOOST_CHECK(2 == count);
  }
}
