/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE simple_pore test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/SimplePore.hpp>
#include <utils/Vector.hpp>

/*
 * The pore has mirror symmetry in z about its center (the documented
 * invariant in SimplePore::calculate_dist). Two probe points at the +z and
 * -z mouths with identical r and |z| therefore MUST yield the same signed
 * distance, and a point inside the solid wall material at either mouth must
 * report a negative (inside) distance.
 *
 * Bug #13: the smoothing-area side test used the signed z in
 * Utils::sqr(z - c_z) instead of Utils::sqr(std::abs(z) - c_z). At the
 * negative-z mouth this forces side=+1 (outside) even for a point inside the
 * solid, flipping the sign of the distance.
 */
BOOST_AUTO_TEST_CASE(mirror_symmetry_in_z) {
  Shapes::SimplePore pore;
  pore.set_axis({0., 0., 1.});
  pore.set_radius(3.);
  pore.set_smoothing_radius(0.1);
  pore.set_length(6.);
  pore.center() = Utils::Vector3d{10., 10., 10.};
  /* derived: c_r = 3.1, c_z = 2.9 */

  /* Probe r = 3.05 (x-offset), z = +/-2.95 relative to the center: inside the
     solid wall material, inside the smoothing fillet at each mouth. r is
     chosen slightly larger than m_rad so the smoothing branch is exercised. */
  Utils::Vector3d const center{10., 10., 10.};
  auto const pos_plus = center + Utils::Vector3d{3.05, 0., 2.95};
  auto const pos_minus = center + Utils::Vector3d{3.05, 0., -2.95};

  double dist_plus, dist_minus;
  Utils::Vector3d vec_plus, vec_minus;
  pore.calculate_dist(pos_plus, dist_plus, vec_plus);
  pore.calculate_dist(pos_minus, dist_minus, vec_minus);

  /* The signed scalar distance must be mirror-symmetric in z. */
  BOOST_TEST_INFO("dist(+z)=" << dist_plus << " dist(-z)=" << dist_minus);
  BOOST_REQUIRE_CLOSE(dist_plus, dist_minus, 1e-9);

  /* The probe lies inside the solid wall material, so both mouths must
     report a negative (inside) distance. */
  BOOST_REQUIRE_LT(dist_plus, 0.);
  BOOST_REQUIRE_LT(dist_minus, 0.);
}
