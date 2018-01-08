/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <limits>

#define BOOST_TEST_MODULE ellipsoid test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "shapes/Ellipsoid.hpp"

bool check_distance_function(const Shapes::Shape &s) {
  double semiaxes[3] = {3.1, 2.2, 1.3};

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      double theta = acos(2. * j * 1e-2 - 1.);

      double dist[3];
      double d;

      /* pos part of surface */
      double pos[3] = {semiaxes[0] * sqrt(1. - i * i * 1e-4) * cos(theta),
                       semiaxes[1] * sqrt(1. - i * i * 1e-4) * sin(theta),
                       semiaxes[2] * i * 1e-2};

      /* check that points on ellipsoid yield zero distance */
      s.calculate_dist(pos, &d, dist);
      if (std::abs(d) > 1e-12)
        return false;

      /* pos outside of surface */
      for (int dim = 0; dim < 3; dim++)
        pos[dim] += 2.3 - dim;
      s.calculate_dist(pos, &d, dist);

      /* trivial test */
      if ((dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2] - d * d) >
          1e-12)
        return false;
    }
  }

  return true;
}

BOOST_AUTO_TEST_CASE(dist_function) {
  Shapes::Ellipsoid e;
  e.semiaxis_a() = 3.1;
  e.semiaxis_b() = 2.2;
  e.semiaxis_c() = 1.3;

  BOOST_CHECK(check_distance_function(e));
}
