/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <shapes/Shape.hpp>
#include <shapes/Wall.hpp>

#include <utils/Vector.hpp>

BOOST_AUTO_TEST_CASE(dist_function) {
  Shapes::Wall shape;
  shape.set_normal(Utils::Vector3d{3., 5., 7.});
  shape.d() = 0.2;

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      for (int k = 0; k < 100; k++) {
        Utils::Vector3d pos = {i * 0.1, j * 0.1, k * 0.1};
        Utils::Vector3d dist{};
        double d;

        shape.calculate_dist(pos, d, dist);
        BOOST_REQUIRE_LE((dist.norm2() - d * d), 1e-12);
        BOOST_REQUIRE_EQUAL(shape.is_inside(pos), d <= 0.);

        /* check if the returned closest point really is on the surface */
        pos -= dist;

        shape.calculate_dist(pos, d, dist);
        BOOST_REQUIRE_LE(std::abs(d), 1e-12);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(rasterize_function) {
  {
    Shapes::Wall shape;
    shape.set_normal(Utils::Vector3d{1., 0., 0.});
    shape.d() = 1.0;
    auto const agrid = 1.0;

    auto const raster = shape.rasterize({5, 5, 5}, agrid, 0.5);
    for (int i = 0; i < 25; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 1);
    }
    for (int i = 25; i < 125; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 0);
    }
  }
  // edge case: wall right before the second slice of LB nodes
  {
    Shapes::Wall shape;
    shape.set_normal(Utils::Vector3d{1., 0., 0.});
    shape.d() = 1.49999999;
    auto const agrid = 1.0;

    auto const raster = shape.rasterize({5, 5, 5}, agrid, 0.5);
    for (int i = 0; i < 25; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 1);
    }
    for (int i = 25; i < 125; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 0);
    }
  }
  // edge case: wall right on the second slice of LB nodes
  {
    Shapes::Wall shape;
    shape.set_normal(Utils::Vector3d{1., 0., 0.});
    shape.d() = 1.50000000;
    auto const agrid = 1.0;

    auto const raster = shape.rasterize({5, 5, 5}, agrid, 0.5);
    for (int i = 0; i < 2 * 25; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 1);
    }
    for (int i = 2 * 25; i < 125; ++i) {
      BOOST_REQUIRE_EQUAL(raster[i], 0);
    }
  }
}
