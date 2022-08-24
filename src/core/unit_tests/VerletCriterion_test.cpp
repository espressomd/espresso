/*
 * Copyright (C) 2022 The ESPResSo project
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
#include "cell_system/CellStructure.hpp"
#include "config/config.hpp"
#include "nonbonded_interactions/VerletCriterion.hpp"
#include "particle_data.hpp"

BOOST_AUTO_TEST_CASE(VerletCriterion_test) {
  auto constexpr skin = 0.4;
  auto constexpr max_cut = 2.5;
  auto constexpr coulomb_cut = 2.0;
  auto constexpr dipolar_cut = 1.8;
  auto constexpr collision_cut = 1.6;

  struct GetMaxCutoff {
    double operator()(int, int) const { return skin + max_cut; }
  };
  struct GetZeroCutoff {
    double operator()(int, int) const { return -skin; }
  };

  VerletCriterion<GetMaxCutoff> criterion(skin, max_cut);
  VerletCriterion<GetMaxCutoff> criterion_inactive(skin, INACTIVE_CUTOFF);
  VerletCriterion<GetZeroCutoff> criterion_long_range(
      skin, max_cut, coulomb_cut, dipolar_cut, collision_cut);

  Particle p1, p2;
  p1.id() = 1;
  p2.id() = 2;

  {
    auto constexpr cutoff = skin + max_cut;
    auto const below = Distance{Utils::Vector3d{cutoff - 0.1, 0.0, 0.0}};
    auto const above = Distance{Utils::Vector3d{cutoff + 0.1, 0.0, 0.0}};
    BOOST_CHECK(criterion(p1, p2, below));
    BOOST_CHECK(!criterion_inactive(p1, p2, below));
    BOOST_CHECK(!criterion(p1, p2, above));
    BOOST_CHECK(!criterion_inactive(p1, p2, above));
  }

#ifdef ELECTROSTATICS
  {
    auto constexpr cutoff = skin + coulomb_cut;
    auto const below = Distance{Utils::Vector3d{cutoff - 0.1, 0.0, 0.0}};
    auto const above = Distance{Utils::Vector3d{cutoff + 0.1, 0.0, 0.0}};
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
#endif // ELECTROSTATICS

#ifdef DIPOLES
  {
    auto constexpr cutoff = skin + dipolar_cut;
    auto const below = Distance{Utils::Vector3d{cutoff - 0.1, 0.0, 0.0}};
    auto const above = Distance{Utils::Vector3d{cutoff + 0.1, 0.0, 0.0}};
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
#endif // DIPOLES

#ifdef COLLISION_DETECTION
  {
    auto constexpr cutoff = skin + collision_cut;
    auto const below = Distance{Utils::Vector3d{cutoff - 0.1, 0.0, 0.0}};
    auto const above = Distance{Utils::Vector3d{cutoff + 0.1, 0.0, 0.0}};
    BOOST_CHECK(criterion_long_range(p1, p2, below));
    BOOST_CHECK(!criterion_long_range(p1, p2, above));
  }
#endif // COLLISION_DETECTION
}
