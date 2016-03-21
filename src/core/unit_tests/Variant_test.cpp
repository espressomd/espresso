/*
  Copyright (C) 2010,2011,2012,2013,2014,2015 The ESPResSo project
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
#include <boost/test/included/unit_test.hpp>

#include "../Variant.hpp"

#include <algorithm>
#include <iostream>

/** Number of nontrivial Baxter permutations of length 2n-1. (A001185) */
#define TEST_NUMBERS { 0, 1, 1, 7, 21, 112, 456, 2603, 13203 }
#define TEST_NUMBERS_PARTIAL_NORM2 { 0, 1, 2, 51, 492, 13036 }
const int test_numbers[] = TEST_NUMBERS;
const int test_numbers_partial_norm2[] = TEST_NUMBERS_PARTIAL_NORM2;
const int n_test_numbers = sizeof(test_numbers) / sizeof(int);


BOOST_AUTO_TEST_CASE(test_constructor) {
  std::vector<int> iv(n_test_numbers), jv;

  std::copy(test_numbers, test_numbers + n_test_numbers, iv.begin());

  Variant v(iv);

  jv = v;

  BOOST_CHECK(iv == jv);
  BOOST_CHECK(v.type() == Variant::INT_VECTOR);

  v = jv;

  BOOST_CHECK(iv == static_cast<std::vector<int> &>(v));
}
