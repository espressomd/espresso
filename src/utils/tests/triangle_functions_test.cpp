/*
 * Copyright (C) 2019 The ESPResSo project
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

/* Unit test for Utils triangle algorithms. */

#define BOOST_TEST_MODULE Utils::triangle_functions
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/triangle_functions.hpp"

auto const epsilon = std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_CASE(normal_) {
  /* Non-degenerate case */
  {
    auto const P1 = Utils::Vector3d{0, 0, 0};
    auto const P2 = Utils::Vector3d{0, 1, 0};
    auto const P3 = Utils::Vector3d{0, 0, 1};

    auto const normal = Utils::get_n_triangle(P1, P2, P3);

    auto const P1P2 = P2 - P1;
    auto const P1P3 = P3 - P1;

    /* Normal should be orthogonal to the sides of the triangle */
    BOOST_CHECK_SMALL(normal * P1P2, epsilon);
    BOOST_CHECK_SMALL(normal * P1P3, epsilon);

    /* Orientation */
    BOOST_CHECK(vector_product(normal, P1P2) * P1P3 > 0.);
    BOOST_CHECK(vector_product(P1P3, normal) * P1P2 > 0.);
    BOOST_CHECK_SMALL(((normal + Utils::get_n_triangle(P2, P1, P3)).norm()),
                      epsilon);
  }

  /* Degenerate case */
  {
    auto const P1 = Utils::Vector3d{-1, 0, 0};
    auto const P2 = Utils::Vector3d{0, 0, 0};
    auto const P3 = Utils::Vector3d{1, 0, 0};

    auto const normal = Utils::get_n_triangle(P1, P2, P3);

    auto const P1P2 = P2 - P1;
    auto const P1P3 = P3 - P1;

    /* Normal should be orthogonal to the sides of the triangle */
    BOOST_CHECK_SMALL(normal * P1P2, epsilon);
    BOOST_CHECK_SMALL(normal * P1P3, epsilon);
  }
}

BOOST_AUTO_TEST_CASE(angle_triangles) {
  /* Notes from @icimrak from Github #3385
  As an example, consider 4 points A,B,C,D in space given by coordinates A =
  [1,1,1], B = [2,1,1], C = [1,2,1], D = [1,1,2]. We want to determine the angle
  between triangles ABC and ACD. In case that the orientations of the triangle
  ABC is [0,0,1] and orientation of ACD is [1,0,0], then the resulting angle
  must be PI/2.0. To get correct result, note that the common edge is AC, and
  one must call the method as angle_btw_triangles(B,A,C,D). With this call we
  have ensured that N1 = AB x AC (which coincides with [0,0,1]) and N2 = AC x AD
  (which coincides with [1,0,0]). Alternatively, if the orientations of the two
  triangles were the opposite, the correct call would be
  angle_btw_triangles(B,C,A,D) so that N1 = CB x CA and N2 = CA x CD.
  */

  const Utils::Vector3d a{1, 1, 1}, b{2, 1, 1}, c{1, 2, 1}, d{1, 1, 2};
  using Utils::angle_btw_triangles;
  BOOST_CHECK_SMALL(std::abs(angle_btw_triangles(b, a, c, d) - M_PI / 2.0),
                    epsilon);
  BOOST_CHECK_SMALL(std::abs(angle_btw_triangles(b, c, a, d) - 3 * M_PI / 2.0),
                    epsilon);
}
