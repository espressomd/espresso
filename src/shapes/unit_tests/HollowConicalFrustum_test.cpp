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

#define BOOST_TEST_MODULE Cone test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/HollowConicalFrustum.hpp>

#include <utils/Vector.hpp>
#include <utils/math/cylindrical_transformation_parameters.hpp>

#include <limits>

BOOST_AUTO_TEST_CASE(dist_function) {
  constexpr double L = 8.0;
  constexpr double R1 = 2.0;
  constexpr double R2 = 3.0;

  constexpr double eps = 10 * std::numeric_limits<double>::epsilon();

  {
    auto ctp = std::make_shared<Utils::CylindricalTransformationParameters>();
    Shapes::HollowConicalFrustum c(R1, R2, L, 0., 1, 0., ctp);

    auto pos = Utils::Vector3d{0.0, 0.0, L / 2.0};
    Utils::Vector3d vec;
    double dist;

    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist - 2.0, eps);
    BOOST_CHECK_SMALL(dist - vec.norm(), eps);

    pos = {{R1, 0.0, L / 2.0}};
    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist, eps);
    BOOST_CHECK_SMALL(dist - vec.norm(), eps);

    pos = {{3.0, 0.0, -L / 2.0}};
    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist, eps);
    BOOST_CHECK_SMALL(dist - vec.norm(), eps);

    c.set_thickness(1.0);
    c.set_r2(R1);
    pos = {{R1 + 1.0, 0.0, L / 2.0}};
    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist - .5, eps);
  }
  {
    auto ctp = std::make_shared<Utils::CylindricalTransformationParameters>(
        Utils::Vector3d{{0., 0., 0.}}, Utils::Vector3d{{1., 0., 0.}});
    Shapes::HollowConicalFrustum c(R1, R2, L, 0., 1, 0., ctp);

    auto pos = Utils::Vector3d{L / 2.0, 0.0, 0.0};
    Utils::Vector3d vec;
    double dist;

    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist - 2.0, eps);
    BOOST_CHECK_SMALL(dist - vec.norm(), eps);

    pos = {{L / 2.0, R1, 0.0}};
    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist, eps);
    BOOST_CHECK_CLOSE(dist, vec.norm(), eps);

    pos = {{-L / 2.0, R2, 0.0}};
    c.calculate_dist(pos, dist, vec);
    BOOST_CHECK_SMALL(dist, eps);
    BOOST_CHECK_SMALL(dist - vec.norm(), eps);
  }
}
