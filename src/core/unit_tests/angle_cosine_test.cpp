/*
 * Copyright (C) 2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "AngleCosine bond"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "bonded_interactions/angle_cosine.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <numbers>

/* Collinear configurations have sin(phi) = 0 in the denominator of the force
 * factor. The forces must be finite and vanish by symmetry. */
BOOST_AUTO_TEST_CASE(angle_cosine_collinear_force_finite) {
  AngleCosineBond const bond(1., std::numbers::pi);

  // phi = pi (stretched) and phi = 0 (folded)
  auto const vec1 = Utils::Vector3d{{-1., 0., 0.}};
  for (auto const &vec2 :
       {Utils::Vector3d{{1., 0., 0.}}, Utils::Vector3d{{-2., 0., 0.}}}) {
    auto const [f_mid, f_left, f_right] = bond.forces(vec1, vec2);
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK(std::isfinite(f_mid[i]));
      BOOST_CHECK(std::isfinite(f_left[i]));
      BOOST_CHECK(std::isfinite(f_right[i]));
      BOOST_CHECK_SMALL(f_mid[i], 1e-9);
      BOOST_CHECK_SMALL(f_left[i], 1e-9);
      BOOST_CHECK_SMALL(f_right[i], 1e-9);
    }
  }
}
