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

#define BOOST_TEST_MODULE "FENE bond"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "bonded_interactions/fene.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <optional>

/* A stretched FENE bond inside the allowed range returns the analytic
 * energy and a non-empty force. This guards the common case and makes
 * sure the compression fix below does not break ordinary evaluation. */
BOOST_AUTO_TEST_CASE(fene_matches_analytic) {
  auto const k = 2.0;
  auto const drmax = 1.5;
  auto const r0 = 2.0;
  FeneBond const fene{k, drmax, r0};

  // stretched bond: len = 3.0 -> dr = +1.0, |dr| < drmax
  auto const dx_stretch = Utils::Vector3d{{3.0, 0.0, 0.0}};
  auto const dr_stretch = dx_stretch.norm() - r0;
  auto const expected_stretch =
      -0.5 * k * drmax * drmax *
      std::log(1.0 - dr_stretch * dr_stretch / (drmax * drmax));
  auto const energy_stretch = fene.energy(dx_stretch);
  BOOST_REQUIRE(energy_stretch.has_value());
  BOOST_CHECK_CLOSE(*energy_stretch, expected_stretch, 1e-10);
  BOOST_CHECK(fene.force(dx_stretch).has_value());

  // compressed bond still within range: len = 1.0 -> dr = -1.0, |dr| < drmax
  auto const dx_compress = Utils::Vector3d{{1.0, 0.0, 0.0}};
  auto const dr_compress = dx_compress.norm() - r0;
  auto const expected_compress =
      -0.5 * k * drmax * drmax *
      std::log(1.0 - dr_compress * dr_compress / (drmax * drmax));
  auto const energy_compress = fene.energy(dx_compress);
  BOOST_REQUIRE(energy_compress.has_value());
  BOOST_CHECK_CLOSE(*energy_compress, expected_compress, 1e-10);
  BOOST_CHECK(fene.force(dx_compress).has_value());
}

/* FENE energy must be symmetric with the force under compression.
 * r_0 > d_r_max is a supported configuration. For a bond compressed
 * below r_0 - d_r_max we have dr <= -drmax, so the bond is broken
 * and both force() and energy() must return an empty optional. */
BOOST_AUTO_TEST_CASE(fene_symmetric_break) {
  auto const k = 1.0;
  auto const drmax = 1.0;
  auto const r0 = 2.0;
  FeneBond const fene{k, drmax, r0};

  // len = 3.5 -> dr = 1.5, |dr| = 1.5 >= drmax (stretched beyond the range)
  auto const dx_stretch = Utils::Vector3d{{3.5, 0.0, 0.0}};

  auto const force_stretch = fene.force(dx_stretch);
  auto const energy_stretch = fene.energy(dx_stretch);
  BOOST_REQUIRE(!force_stretch.has_value());
  BOOST_REQUIRE(!energy_stretch.has_value());

  // len = 0.5 -> dr = -1.5, |dr| = 1.5 >= drmax (compressed beyond the range)
  auto const dx_compress = Utils::Vector3d{{0.5, 0.0, 0.0}};

  auto const force_compress = fene.force(dx_compress);
  auto const energy_compress = fene.energy(dx_compress);
  BOOST_REQUIRE(!force_compress.has_value());
  BOOST_REQUIRE(!energy_compress.has_value());
}
