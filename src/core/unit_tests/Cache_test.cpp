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

#include "utils/Cache.hpp"
using Utils::Cache;

BOOST_AUTO_TEST_CASE(types) {
  using cache_type = Cache<int, char>;
  static_assert(std::is_same<cache_type::key_type, int>::value, "");
  static_assert(std::is_same<cache_type::value_type, const char *>::value, "");
}

BOOST_AUTO_TEST_CASE(get_value) {
  {
    Cache<int, char> cache;
    cache.put(41, 42);
    cache.put(42, 43);
    BOOST_REQUIRE(cache.get(41));
    BOOST_CHECK(42 == *cache.get(41));
    BOOST_REQUIRE(cache.get(42));
    BOOST_CHECK(43 == *cache.get(42));
  }

  {
    Cache<int, char> cache;
    BOOST_CHECK(nullptr == cache.get(123));
    BOOST_CHECK(nullptr == cache.get(11));
  }
}

BOOST_AUTO_TEST_CASE(has) {
  {
    Cache<int, char> cache;
    cache.put(41,0);
    cache.put(1, 0);

    BOOST_CHECK(cache.has(41));
    BOOST_CHECK(cache.has(1));
    BOOST_CHECK(!cache.has(2));
    BOOST_CHECK(cache.has(1));
  }
}

BOOST_AUTO_TEST_CASE(invalidate) {
  {
    Cache<int, char> cache;
    cache.put(41, 11);
    cache.put(1, 123);

    cache.invalidate();

    BOOST_CHECK(!cache.has(41));
    BOOST_CHECK(!cache.has(1));
    BOOST_CHECK(!cache.has(2));
  }
}

BOOST_AUTO_TEST_CASE(put) {
  Cache<int, char> cache;

  cache.put(0, 2);
}
