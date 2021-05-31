/*
 * Copyright (C) 2020 The ESPResSo project
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

#define BOOST_TEST_MODULE unordered_map test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/serialization/unordered_map.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>
#include <unordered_map>

BOOST_AUTO_TEST_CASE(serialization) {

  std::unordered_map<int, char> const map_original{{65, 'A'}, {66, 'B'}};

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << map_original;

  std::unordered_map<int, char> map_deserialized;
  boost::archive::text_iarchive in_ar(stream);
  in_ar >> map_deserialized;

  BOOST_CHECK_EQUAL(map_deserialized.size(), 2);
  BOOST_CHECK_EQUAL(map_deserialized.count(65), 1);
  BOOST_CHECK_EQUAL(map_deserialized.count(66), 1);
  BOOST_CHECK_EQUAL(map_deserialized.at(65), 'A');
  BOOST_CHECK_EQUAL(map_deserialized.at(66), 'B');
}
