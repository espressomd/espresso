/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#define BOOST_TEST_MODULE acumulator test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Accumulator.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <algorithm>
#include <limits>
#include <sstream>
#include <vector>

BOOST_AUTO_TEST_CASE(accumulator) {
  auto const test_data1 = std::vector<double>{{0.0, 1.0, 2.0, 3.0}};
  auto const test_data2 = std::vector<double>{{1.5, 3.5, 3.5, 4.5}};
  auto const test_mean = std::vector<double>{{0.75, 2.25, 2.75, 3.75}};
  auto const test_var = std::vector<double>{{1.125, 3.125, 1.125, 1.125}};
  auto const test_stderr = std::vector<double>{{0.75, 1.25, 0.75, 0.75}};

  // check statistics
  auto acc = Utils::Accumulator(4);
  acc(test_data1);
  BOOST_CHECK(acc.mean() == test_data1);
  BOOST_CHECK(acc.variance() ==
              std::vector<double>(4, std::numeric_limits<double>::max()));
  acc(test_data2);
  BOOST_CHECK((acc.mean() == test_mean));
  BOOST_CHECK((acc.variance() == test_var));
  BOOST_CHECK((acc.std_error() == test_stderr));
  BOOST_CHECK_THROW(acc(std::vector<double>{{1, 2, 3, 4, 5}}),
                    std::runtime_error);

  // check serialization
  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << acc;

  auto acc_deserialized = Utils::Accumulator(4);
  boost::archive::text_iarchive in_ar(stream);
  in_ar >> acc_deserialized;

  BOOST_CHECK((acc_deserialized.mean() == test_mean));
  BOOST_CHECK((acc_deserialized.variance() == test_var));
  BOOST_CHECK((acc_deserialized.std_error() == test_stderr));
}
