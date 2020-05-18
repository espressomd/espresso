/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Wall test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/Wall.hpp>

bool check_distance_function(const Shapes::Shape &s) {
  for (int i = 0; i < 100; i++)
    for (int j = 0; j < 100; j++)
      for (int k = 0; k < 100; k++) {
        Utils::Vector3d pos = {i * 0.1, j * 0.1, k * 0.1};
        Utils::Vector3d dist{};
        double d;

        s.calculate_dist(pos, d, dist);
        /* trivial test */
        if ((dist.norm2() - d * d) > 1e-12) {
          return false;
        }

        /* check if the returned closest point really is on the surface */
        pos -= dist;

        s.calculate_dist(pos, d, dist);
        if (std::abs(d) > 1e-12) {
          return false;
        }
      }

  return true;
}

BOOST_AUTO_TEST_CASE(dist_function) {
  Shapes::Wall w;
  w.set_normal(Utils::Vector3d{3., 5., 7.});
  w.d() = 0.2;

  BOOST_CHECK(check_distance_function(w));
}

/* @todo  Functional unit test of the distance function */
