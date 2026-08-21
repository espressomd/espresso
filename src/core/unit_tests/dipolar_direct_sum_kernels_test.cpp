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

#define BOOST_TEST_MODULE "dipolar direct sum kernels"
#define BOOST_TEST_DYN_LINK

#include <config/config.hpp>

#ifdef ESPRESSO_DIPOLES

#include <boost/test/unit_test.hpp>

#include "magnetostatics/dipolar_direct_sum_kernels.hpp"

#include <utils/Vector.hpp>

BOOST_AUTO_TEST_SUITE(suite)

// Anchor: two identical dipoles m=(1,0,0) separated by d=(1,0,0).
// Analytically f=(-6,0,0), torque=0 (parallel moments => no torque).
BOOST_AUTO_TEST_CASE(pair_force_analytic) {
  Utils::Vector3d const d{1., 0., 0.};
  Utils::Vector3d const m{1., 0., 0.};
  auto const pf = pair_force(d, m, m);
  BOOST_CHECK_CLOSE(pf.f[0], -6., 1e-10);
  BOOST_CHECK_SMALL(pf.f[1], 1e-12);
  BOOST_CHECK_SMALL(pf.f[2], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[0], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[1], 1e-12);
  BOOST_CHECK_SMALL(pf.torque[2], 1e-12);
}

// Anchor: two identical dipoles m=(1,0,0) separated by d=(1,0,0).
// U = (m1.m2)/r^3 - 3 (m1.d)(m2.d)/r^5 = 1 - 3 = -2.
BOOST_AUTO_TEST_CASE(pair_potential_analytic) {
  Utils::Vector3d const d{1., 0., 0.};
  Utils::Vector3d const m{1., 0., 0.};
  BOOST_CHECK_CLOSE(pair_potential(d, m, m), -2., 1e-10);

  // Perpendicular moment: m2.d = 0, so U = (m1.m2)/r^3 = 0.
  Utils::Vector3d const mp{0., 1., 0.};
  BOOST_CHECK_SMALL(pair_potential(d, m, mp), 1e-12);
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
// Anchor: field of dipole m=(1,0,0) at d=(1,0,0).
// H = 3 (m.d) d / r^5 - m / r^3 = 3*(1,0,0) - (1,0,0) = (2,0,0).
BOOST_AUTO_TEST_CASE(dipole_field_analytic) {
  Utils::Vector3d const d{1., 0., 0.};
  Utils::Vector3d const m{1., 0., 0.};
  auto const h = dipole_field(d, m);
  BOOST_CHECK_CLOSE(h[0], 2., 1e-10);
  BOOST_CHECK_SMALL(h[1], 1e-12);
  BOOST_CHECK_SMALL(h[2], 1e-12);
}
#endif // ESPRESSO_DIPOLE_FIELD_TRACKING

BOOST_AUTO_TEST_SUITE_END()

#else  // ESPRESSO_DIPOLES
int main(int argc, char **argv) {}
#endif // ESPRESSO_DIPOLES
