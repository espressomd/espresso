/*
Copyright (C) 2010-2018 The ESPResSo project

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
