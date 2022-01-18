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

#define BOOST_TEST_MODULE NoWhere test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/NoWhere.hpp>
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <limits>

BOOST_AUTO_TEST_CASE(dist_function) {
  constexpr auto infinity = std::numeric_limits<double>::infinity();

  Shapes::NoWhere const s;

  Utils::Vector3d const positions[2] = {
      {0.0, 1.0, 2.0},
      {-10.0, 0.1, 5.0},
  };

  for (auto const &pos : positions) {
    Utils::Vector3d dist{};
    double d;

    s.calculate_dist(pos, d, dist);
    BOOST_CHECK_EQUAL(d, infinity);

    for (auto xyz : dist) {
      BOOST_CHECK_EQUAL(xyz, infinity);
    }
  }
}
