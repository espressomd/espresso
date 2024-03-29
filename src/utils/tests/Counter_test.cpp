/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE Utils::Counter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Counter.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>

using Utils::Counter;

BOOST_AUTO_TEST_CASE(ctor) {
  {
    auto c = Counter<int>(5);
    BOOST_CHECK_EQUAL(c.value(), 5);
    BOOST_CHECK_EQUAL(c.initial_value(), 5);
  }

  {
    auto c = Counter<int>(5, 6);
    BOOST_CHECK_EQUAL(c.initial_value(), 5);
    BOOST_CHECK_EQUAL(c.value(), 6);
  }

  {
    auto c = Counter<int>();
    BOOST_CHECK_EQUAL(c.initial_value(), 0);
    BOOST_CHECK_EQUAL(c.value(), 0);
  }
}

BOOST_AUTO_TEST_CASE(increment) {
  auto c = Counter<int>(5);
  c.increment();
  BOOST_CHECK_EQUAL(c.value(), 6);
  BOOST_CHECK_EQUAL(c.initial_value(), 5);
}

BOOST_AUTO_TEST_CASE(serialization) {
  auto c = Counter<int>(5);
  c.increment();

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << c;

  boost::archive::text_iarchive in_ar(stream);
  Counter<int> c_deserialized;
  in_ar >> c_deserialized;

  BOOST_CHECK_EQUAL(c_deserialized.value(), 6);
  BOOST_CHECK_EQUAL(c_deserialized.initial_value(), 5);
}
