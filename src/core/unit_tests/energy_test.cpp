/*
 * Copyright (C) 2021 The ESPResSo project
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

#define BOOST_TEST_MODULE tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "energy_inline.hpp"

BOOST_AUTO_TEST_CASE(translational_kinetic_energy_) {

  Particle p;
#ifdef MASS
  p.p.mass = 2.;
#endif
  p.m.v = {3., 4., 5.};

  auto const expected = 0.5 * p.p.mass * p.m.v.norm2();
  BOOST_TEST(translational_kinetic_energy(p) == expected);

/* virtual */
#ifdef VIRTUAL_SITES
  {
    Particle p;
#ifdef MASS
    p.p.mass = 2.;
#endif
    p.p.is_virtual = true;
    p.m.v = {3., 4., 5.};

    auto const expected = 0.;
    BOOST_TEST(translational_kinetic_energy(p) == expected);
  }
#endif
}

BOOST_AUTO_TEST_CASE(rotational_kinetic_energy_) {
  BOOST_TEST(rotational_kinetic_energy(Particle{}) == 0.);

#ifdef ROTATION
  {
    Particle p;
    p.m.omega = {1., 2., 3.};
    p.p.rotation = 1;

    auto const expected =
        0.5 * (hadamard_product(p.m.omega, p.m.omega) * p.p.rinertia);
    BOOST_TEST(rotational_kinetic_energy(p) == expected);
  }
#endif
}

BOOST_AUTO_TEST_CASE(kinetic_energy_) {
  Particle p;
#ifdef MASS
  p.p.mass = 2.;
#endif
  p.m.v = {3., 4., 5.};

#ifdef ROTATION
  p.m.omega = {1., 2., 3.};
  p.p.rotation = 1;
#endif

  auto const expected =
      translational_kinetic_energy(p) + rotational_kinetic_energy(p);
  BOOST_TEST(calc_kinetic_energy(p) == expected);
}