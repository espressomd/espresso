/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE Verlet criterion checks
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "ParticleStoreTestFixture.hpp"
#include "config/config.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "system/System.hpp"

#include <utils/math/sqr.hpp>

BOOST_AUTO_TEST_CASE(VerletCriterion_test) {
  auto constexpr skin = 0.4;
  auto constexpr max_cut = 2.5;
  auto constexpr coulomb_cut = 2.0;
  auto constexpr dipolar_cut = 1.8;
  auto constexpr collision_cut = 1.6;

  struct GetMaxCutoff {
    GetMaxCutoff(System::System const &) {}
    double operator()(int, int) const { return skin + max_cut; }
  };
  struct GetZeroCutoff {
    GetZeroCutoff(System::System const &) {}
    double operator()(int, int) const { return -skin; }
  };

  auto const &system = System::get_system();
  VerletCriterion<GetMaxCutoff> criterion(system, skin, max_cut);
  VerletCriterion<GetMaxCutoff> criterion_inactive(system, skin,
                                                   inactive_cutoff);
  VerletCriterion<GetZeroCutoff> criterion_long_range(
      system, skin, max_cut, coulomb_cut, dipolar_cut, collision_cut);

  // id/charge/dipole moment live in the ParticleStore; attach both hand-made
  // particles to a standalone store. The VerletCriterion only reads these
  // fields by const ref, so a standalone fixture is sufficient.
  ParticleStoreTestFixture fixture{};
  auto p1 = fixture.make();
  auto p2 = fixture.make();
  p1.id() = 1;
  p2.id() = 2;

  {
    auto constexpr cutoff = skin + max_cut;
    auto const below = Utils::sqr(cutoff - 0.1);
    auto const above = Utils::sqr(cutoff + 0.1);
    BOOST_CHECK(criterion(p1, p2, below));
    BOOST_CHECK(!criterion_inactive(p1, p2, below));
    BOOST_CHECK(!criterion(p1, p2, above));
    BOOST_CHECK(!criterion_inactive(p1, p2, above));
  }

#ifdef ESPRESSO_ELECTROSTATICS
  {
    auto constexpr cutoff = skin + coulomb_cut;
    auto const below = Utils::sqr(cutoff - 0.1);
    auto const above = Utils::sqr(cutoff + 0.1);
    BOOST_CHECK(!criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p2.q() = 1.;
    BOOST_CHECK(!criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p1.q() = 1.;
    BOOST_CHECK(criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p1.q() = 0.;
    p2.q() = 0.;
  }
#endif // ESPRESSO_ELECTROSTATICS

#ifdef ESPRESSO_DIPOLES
  {
    auto constexpr cutoff = skin + dipolar_cut;
    auto const below = Utils::sqr(cutoff - 0.1);
    auto const above = Utils::sqr(cutoff + 0.1);
    BOOST_CHECK(!criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p2.dipm() = 1.;
    BOOST_CHECK(!criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p1.dipm() = 1.;
    BOOST_CHECK(criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
    p1.dipm() = 0.;
    p2.dipm() = 0.;
  }
#endif // ESPRESSO_DIPOLES

#ifdef ESPRESSO_COLLISION_DETECTION
  {
    auto constexpr cutoff = skin + collision_cut;
    auto const below = Utils::sqr(cutoff - 0.1);
    auto const above = Utils::sqr(cutoff + 0.1);
    BOOST_CHECK(criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
  }
#endif // ESPRESSO_COLLISION_DETECTION
}
