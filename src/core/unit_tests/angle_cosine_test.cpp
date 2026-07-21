/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE angle_cosine
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "bonded_interactions/angle_cosine.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <tuple>

/** Regression test: at the collinear equilibrium (phi0=pi), forces must be
 *  finite (not NaN / inf).  Before the fix, angle_generic_force was called
 *  with sanitize_cosine=false, so cos_phi could reach exactly -1, making
 *  sin_phi = sqrt(1 - cos²) = 0 and the force factor diverge.
 */
BOOST_AUTO_TEST_CASE(angle_cosine_collinear_force_finite) {
  // bend=1.0, phi0=pi (the default equilibrium — collinear configuration)
  AngleCosineBond bond(1.0, M_PI);

  // Three collinear particles: p1-p2-p3 along the x axis.
  // vec1 = p1 - p2 = {-1, 0, 0}  (central particle is p2)
  // vec2 = p3 - p2 = { 1, 0, 0}
  Utils::Vector3d const vec1 = {-1., 0., 0.};
  Utils::Vector3d const vec2 = {1., 0., 0.};

  auto const [f_mid, f_left, f_right] = bond.forces(vec1, vec2);

  // All force components must be finite (no NaN, no inf).
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK(std::isfinite(f_mid[i]));
    BOOST_CHECK(std::isfinite(f_left[i]));
    BOOST_CHECK(std::isfinite(f_right[i]));
  }

  // At exact equilibrium the net force on every particle is zero (or very
  // small — the sanitized cosine is clamped to ~0.9999999999, not exactly -1,
  // so a tiny residual of order 1e-10 is acceptable).
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(f_mid[i], 1e-9);
    BOOST_CHECK_SMALL(f_left[i], 1e-9);
    BOOST_CHECK_SMALL(f_right[i], 1e-9);
  }
}
