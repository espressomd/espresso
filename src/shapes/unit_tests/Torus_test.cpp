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

#define BOOST_TEST_MODULE torus test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/Torus.hpp>
#include <utils/Vector.hpp>

#include <cmath>

BOOST_AUTO_TEST_CASE(torus_axis_no_nan) {
  Shapes::Torus t;
  t.m_center = {5., 5., 5.};
  t.set_normal({0., 0., 1.});
  t.set_radius(3.);
  t.set_tube_radius(1.);
  t.direction() = 1.0;

  /* Point on the torus symmetry axis: r == 0 in torus frame */
  Utils::Vector3d const pos{5., 5., 4.};
  double dist{};
  Utils::Vector3d vec{};

  t.calculate_dist(pos, dist, vec);

  BOOST_CHECK(std::isfinite(dist));
  BOOST_CHECK(std::isfinite(vec[0]));
  BOOST_CHECK(std::isfinite(vec[1]));
  BOOST_CHECK(std::isfinite(vec[2]));
}
