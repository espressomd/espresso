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

#define BOOST_TEST_MODULE inner_product test
#define BOOST_TEST_DYN_LINK
#include <boost/range/numeric.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/inner_product.hpp"

BOOST_AUTO_TEST_CASE(inner_product) {
  const std::array<int, 3> left_array{1, 2, 9};
  const std::array<double, 3> right_array{0.5, 1.25, 3.1};
  BOOST_CHECK(Utils::inner_product(left_array, right_array) ==
              boost::inner_product(left_array, right_array, 0.0));
}
