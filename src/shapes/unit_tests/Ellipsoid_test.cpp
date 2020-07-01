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

#include <cmath>

#define BOOST_TEST_MODULE ellipsoid test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/Ellipsoid.hpp>

bool check_distance_function(const Shapes::Shape &s) {
  double semiaxes[3] = {3.1, 2.2, 1.3};

  int N = 100;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double theta = 2. * i / N * M_PI;
      double v = j / (N - 1.);

      Utils::Vector3d dist;
      double d;

      /* pos part of surface */
      Utils::Vector3d pos = {semiaxes[0] * sqrt(1. - v * v) * cos(theta),
                             semiaxes[1] * sqrt(1. - v * v) * sin(theta),
                             semiaxes[2] * v};

      /* check that points on ellipsoid yield zero distance */
      s.calculate_dist(pos, d, dist);
      if (std::abs(d) > 1e-12)
        return false;

      /* pos outside of surface */
      for (int dim = 0; dim < 3; dim++)
        pos[dim] += 2.3 - dim;
      s.calculate_dist(pos, d, dist);

      /* trivial test */
      if ((dist.norm2() - d * d) > 1e-12)
        return false;
    }
  }

  return true;
}

BOOST_AUTO_TEST_CASE(dist_function) {
  Shapes::Ellipsoid e;
  e.set_semiaxis_a(3.1);
  e.set_semiaxis_b(2.2);
  e.set_semiaxis_c(1.3);

  BOOST_CHECK(check_distance_function(e));
}
