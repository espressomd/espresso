/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include <boost/test/tools/old/interface.hpp>
#define BOOST_TEST_MODULE sphere test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/Sphere.hpp>
#include <utils/Vector.hpp>

#include <limits>

void check_distance_function(Shapes::Sphere &s) {
  Utils::Vector3d pos;
  Utils::Vector3d vec;
  double dist;
  // multiply by 100 because BOOST_REQUIRE_CLOSE takes a percentage tolerance
  auto const tol = std::numeric_limits<double>::epsilon() * 100;

  s.rad() = 1.0;

  pos = {0., 0., 0.};
  s.calculate_dist(pos, dist, vec);
  double always_pos_dist = -s.direction() * dist;
  BOOST_REQUIRE_GE(always_pos_dist, 0.0);
  BOOST_REQUIRE_CLOSE(always_pos_dist, s.rad(), tol);
  BOOST_REQUIRE_CLOSE(always_pos_dist, vec.norm(), tol);

  for (int i = 0; i < 3; ++i) {
    pos[i] = 1.0;
    s.calculate_dist(pos, dist, vec);
    double always_pos_dist = -s.direction() * dist;
    BOOST_REQUIRE_GE(always_pos_dist, 0.0);
    BOOST_REQUIRE_CLOSE(dist, 0.0, tol);
    BOOST_REQUIRE_CLOSE(always_pos_dist, vec.norm(), tol);
    pos = {0., 0., 0.};
  }
}

BOOST_AUTO_TEST_CASE(dist_function) {
  Shapes::Sphere s_pos;
  Shapes::Sphere s_neg;
  s_pos.direction() = 1.0;
  s_neg.direction() = -1.0;

  check_distance_function(s_pos);
  check_distance_function(s_neg);
}
