/*
 * Copyright (C) 2021-2026 The ESPResSo project
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
#include "ParticleStoreTestFixture.hpp"
#include "PropagationMode.hpp"
#include "energy_inline.hpp"

#include <utils/Vector.hpp>

BOOST_AUTO_TEST_CASE(translational_kinetic_energy_) {
  // real particle
  {
    // Kinematic fields live in the ParticleStore; attach the hand-made
    // particle to a standalone store before setting/reading them.
    ParticleStoreTestFixture fx;
    auto p = fx.make();
#ifdef ESPRESSO_MASS
    p.mass() = 2.;
#endif
    p.v() = {3., 4., 5.};

    auto const expected = 0.5 * p.mass() * p.v().norm2();
    BOOST_CHECK_EQUAL(translational_kinetic_energy(p), expected);
  }

  // virtual particle
  {
#ifdef ESPRESSO_VIRTUAL_SITES

    ParticleStoreTestFixture fx;
    auto p = fx.make();
#ifdef ESPRESSO_MASS
    p.mass() = 2.;
#endif
    p.propagation() = PropagationMode::TRANS_VS_RELATIVE;
    p.v() = {3., 4., 5.};

    auto const expected = 0.;
    BOOST_CHECK_EQUAL(translational_kinetic_energy(p), expected);
#endif
  }
}

BOOST_AUTO_TEST_CASE(rotational_kinetic_energy_) {
  // a default-seeded particle has zero rotation, so its rotational energy is 0
  {
    ParticleStoreTestFixture fx;
    BOOST_CHECK_EQUAL(rotational_kinetic_energy(fx.make()), 0.);
  }

#ifdef ESPRESSO_ROTATION
  {
    ParticleStoreTestFixture fx;
    auto p = fx.make();
    p.omega() = {1., 2., 3.};
    p.set_can_rotate_all_axes();

    auto const expected = 0.5 * (hadamard_product(Utils::Vector3d(p.omega()),
                                                  Utils::Vector3d(p.omega())) *
                                 Utils::Vector3d(p.rinertia()));
    BOOST_CHECK_EQUAL(rotational_kinetic_energy(p), expected);
  }

  // virtual particle
  {
#ifdef ESPRESSO_VIRTUAL_SITES

    ParticleStoreTestFixture fx;
    auto p = fx.make();
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    p.rinertia() = {1., 2., 3.};
#endif
    p.propagation() = PropagationMode::ROT_VS_RELATIVE;
    p.omega() = {3., 4., 5.};
    p.set_can_rotate_all_axes();

    auto const expected = 0.;
    BOOST_CHECK_EQUAL(rotational_kinetic_energy(p), expected);
#endif
  }
#endif
}
