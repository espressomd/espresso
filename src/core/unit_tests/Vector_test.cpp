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

/** \file Vector_test.cpp Unit tests for the Utils::Vector class.
 *
*/

#define BOOST_TEST_MODULE Vector test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <numeric>

#include "../Vector.hpp"

/** Number of nontrivial Baxter permutations of length 2n-1. (A001185) */
#define TEST_NUMBERS                                                           \
  { 0, 1, 1, 7, 21, 112, 456, 2603, 13203 }
#define TEST_NUMBERS_PARTIAL_NORM2                                             \
  { 0, 1, 2, 51, 492, 13036 }
constexpr int test_numbers[] = TEST_NUMBERS;
constexpr int n_test_numbers = sizeof(test_numbers) / sizeof(int);

template <int n> bool il_constructor() {
  bool pass = true;
  Vector<n, int> v(TEST_NUMBERS);

  pass &= v.size() == n;

  for (int i = 0; i < std::min(n_test_numbers, n); i++)
    pass &= v[i] == test_numbers[i];

  return pass;
}

template <int n> bool default_constructor() {
  bool pass = true;
  Vector<n, int> v;

  for (int i = 0; i < n; i++) {
    v[i] = i;
  }

  return pass;
}

template <int n> bool norm2() {
  Vector<n, int> v(std::begin(test_numbers), test_numbers + n);

  return v.norm2() == std::inner_product(v.begin(), v.end(), v.begin(), 0);
}

BOOST_AUTO_TEST_CASE(initializer_list_constructor) {
  Vector<n_test_numbers, int> v(TEST_NUMBERS);

  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));
}

BOOST_AUTO_TEST_CASE(iterator_constructor) {
  Vector<n_test_numbers, int> v(std::begin(test_numbers), std::end(test_numbers));
  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));
}

BOOST_AUTO_TEST_CASE(default_constructor_test) {
  BOOST_CHECK(default_constructor<1>());
  BOOST_CHECK(default_constructor<2>());
  BOOST_CHECK(default_constructor<3>());
  BOOST_CHECK(default_constructor<4>());
  BOOST_CHECK(default_constructor<5>());
  BOOST_CHECK(default_constructor<6>());
  BOOST_CHECK(default_constructor<7>());
  BOOST_CHECK(default_constructor<8>());
  BOOST_CHECK(default_constructor<9>());
  BOOST_CHECK(default_constructor<10>());
}

BOOST_AUTO_TEST_CASE(test_norm2) {
  BOOST_CHECK(norm2<1>());
  BOOST_CHECK(norm2<2>());
  BOOST_CHECK(norm2<3>());
  BOOST_CHECK(norm2<4>());
}
