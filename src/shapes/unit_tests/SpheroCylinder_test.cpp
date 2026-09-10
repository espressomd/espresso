/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#define BOOST_TEST_MODULE spherocylinder test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/SpheroCylinder.hpp>
#include <utils/Vector.hpp>

#include <cmath>

/* Verify that a default-constructed SpheroCylinder has m_direction initialized
 * to 1.0 (outside convention), so calculate_dist returns finite, correctly
 * signed distances without reading garbage memory. */
BOOST_AUTO_TEST_CASE(direction_default_initialized) {
  Shapes::SpheroCylinder sc;
  // Configure geometry without touching direction — that is exactly the
  // scenario that triggered the bug.
  sc.set_radius(1.0);
  sc.set_length(2.0);

  Utils::Vector3d vec;
  double dist;

  // Point clearly outside (radial distance > radius, axial position = 0)
  Utils::Vector3d outside_point{3.0, 0.0, 0.0};
  sc.calculate_dist(outside_point, dist, vec);

  // Distance must be a finite number — not NaN or inf from garbage m_direction
  BOOST_CHECK(std::isfinite(dist));

  // With m_direction == 1.0 (outside convention) a point outside the shape
  // must have positive distance.
  BOOST_CHECK_GT(dist, 0.0);

  // Point clearly inside
  Utils::Vector3d inside_point{0.0, 0.0, 0.0};
  sc.calculate_dist(inside_point, dist, vec);

  BOOST_CHECK(std::isfinite(dist));

  // With m_direction == 1.0 a point inside has negative distance.
  BOOST_CHECK_LT(dist, 0.0);
}
