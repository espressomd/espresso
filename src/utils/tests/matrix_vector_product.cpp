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

#define BOOST_TEST_MODULE matrix_vector_product test
#define BOOST_TEST_DYN_LINK
#include <boost/range/numeric.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/math/matrix_vector_product.hpp"

#include <limits>

extern constexpr std::array<std::array<int, 3>, 3> matrix{
    {{{1, 2, 9}}, {{8, 41, 6}}, {{31, 15, 99}}}};

BOOST_AUTO_TEST_CASE(inner_product) {
  const std::array<double, 3> vector{{0.5, 1.25, 3.1}};
  auto const result = Utils::matrix_vector_product<double, 3, matrix>(vector);
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(result[i], boost::inner_product(matrix[i], vector, 0.0),
                      100. * std::numeric_limits<double>::epsilon());
  }
}
