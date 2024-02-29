/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE energy calculation
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "energy_inline.hpp"

#include <utils/Vector.hpp>

BOOST_AUTO_TEST_CASE(translational_kinetic_energy_) {
  // real particle
  {
    Particle p;
#ifdef MASS
    p.mass() = 2.;
#endif
    p.v() = {3., 4., 5.};

    auto const expected = 0.5 * p.mass() * p.v().norm2();
    BOOST_CHECK_EQUAL(translational_kinetic_energy(p), expected);
  }

  // virtual particle
  {
#ifdef VIRTUAL_SITES

    Particle p;
#ifdef MASS
    p.mass() = 2.;
#endif
    p.set_virtual(true);
    p.v() = {3., 4., 5.};

    auto const expected = 0.;
    BOOST_CHECK_EQUAL(translational_kinetic_energy(p), expected);
#endif
  }
}

BOOST_AUTO_TEST_CASE(rotational_kinetic_energy_) {
  BOOST_CHECK_EQUAL(rotational_kinetic_energy(Particle{}), 0.);

#ifdef ROTATION
  {
    Particle p;
    p.omega() = {1., 2., 3.};
    p.set_can_rotate_all_axes();

    auto const expected =
        0.5 * (hadamard_product(p.omega(), p.omega()) * p.rinertia());
    BOOST_CHECK_EQUAL(rotational_kinetic_energy(p), expected);
  }

  // virtual particle
  {
#ifdef VIRTUAL_SITES

    Particle p;
#ifdef ROTATIONAL_INERTIA
    p.rinertia() = {1., 2., 3.};
#endif
    p.set_virtual(true);
    p.omega() = {3., 4., 5.};
    p.set_can_rotate_all_axes();

    auto const expected = 0.;
    BOOST_CHECK_EQUAL(rotational_kinetic_energy(p), expected);
#endif
  }
#endif
}

BOOST_AUTO_TEST_CASE(kinetic_energy_) {
  Particle p;
#ifdef MASS
  p.mass() = 2.;
#endif
  p.v() = {3., 4., 5.};

#ifdef ROTATION
  p.omega() = {1., 2., 3.};
  p.set_can_rotate_all_axes();
#endif

  auto const expected =
      translational_kinetic_energy(p) + rotational_kinetic_energy(p);
  BOOST_CHECK_EQUAL(calc_kinetic_energy(p), expected);
}
