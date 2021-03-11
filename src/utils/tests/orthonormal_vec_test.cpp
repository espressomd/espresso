/*
 * Copyright (C) 2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Utils::orthonormal_vec test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/math/orthonormal_vec.hpp>

BOOST_AUTO_TEST_CASE(orthonormal_vec_test) {
  constexpr auto eps = 1e-14;

  auto const v0 = Utils::Vector3d{{1.1, -2.2, 3.3}};
  auto v0_orth = Utils::calc_orthonormal_vector(v0);
  BOOST_CHECK_SMALL(v0 * v0_orth, eps);
  BOOST_CHECK_SMALL(1 - v0_orth.norm(), eps);

  auto const v1 = Utils::VectorXd<2>{{1., 0.}};
  auto v1_orth = Utils::calc_orthonormal_vector(v1);
  BOOST_CHECK_SMALL(v1 * v1_orth, eps);
  BOOST_CHECK_SMALL(1 - v1_orth.norm(), eps);
}