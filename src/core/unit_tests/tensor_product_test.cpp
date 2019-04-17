/*
  Copyright (C) 2018 The ESPResSo project

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

/** \file
 * Unit tests for the Utils::tensor_product class.
 *
 */

#define BOOST_TEST_MODULE Utils::tensor_product
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/tensor_product.hpp"
using Utils::tensor_product;

#include <type_traits>

BOOST_AUTO_TEST_CASE(square) {
  Utils::Vector<int, 4> x{1, 2, 3, 4};
  Utils::Vector<int, 4> y{5, 6, 7, 8};

  using expected_type = Utils::Vector<Utils::Vector<int, 4>, 4>;
  using actual_type = decltype(tensor_product(x, y));

  /* Check return type */
  static_assert(std::is_same<expected_type, actual_type>::value, "");

  auto prod = tensor_product(x, y);

  /* Check values */
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::size_t j = 0; j < y.size(); j++) {
      BOOST_CHECK(prod[i][j] == x[i] * y[j]);
    }
}

BOOST_AUTO_TEST_CASE(non_square) {
  Utils::Vector3i x{1, 2, 3};
  Utils::Vector<int, 4> y{5, 6, 7, 8};

  using expected_type = Utils::Vector<Utils::Vector<int, 4>, 3>;
  using actual_type = decltype(tensor_product(x, y));

  /* Check return type */
  static_assert(std::is_same<expected_type, actual_type>::value, "");

  auto prod = tensor_product(x, y);

  /* Check values */
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::size_t j = 0; j < y.size(); j++) {
      BOOST_CHECK(prod[i][j] == x[i] * y[j]);
    }
}

BOOST_AUTO_TEST_CASE(left_scalar) {
  double x = 3.;
  Utils::Vector2d y{1., 2.};

  using expected_type = Utils::Vector2d;
  using actual_type = decltype(tensor_product(x, y));

  /* Check return type */
  static_assert(std::is_same<expected_type, actual_type>::value, "");

  auto prod = tensor_product(x, y);

  BOOST_CHECK((x * y) == prod);
}
