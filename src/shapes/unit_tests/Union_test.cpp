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

#include <memory>

#define BOOST_TEST_MODULE Union test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/Union.hpp>
#include <shapes/Wall.hpp>

BOOST_AUTO_TEST_CASE(dist_function) {
  {
    auto wall1 = std::make_shared<Shapes::Wall>();
    wall1->set_normal(Utils::Vector3d{0., 0., 1.});
    wall1->d() = 0.0;

    auto wall2 = std::make_shared<Shapes::Wall>();
    wall2->set_normal(Utils::Vector3d{0., 0., -1.});
    wall2->d() = -10.0;

    Shapes::Union uni;
    uni.add(wall1);
    uni.add(wall2);

    auto check_union = [&wall1, &wall2, &uni](Utils::Vector3d const &pos) {
      double wall1_dist;
      Utils::Vector3d wall1_vec;
      wall1->calculate_dist(pos, wall1_dist, wall1_vec);

      double wall2_dist;
      Utils::Vector3d wall2_vec;
      wall2->calculate_dist(pos, wall2_dist, wall2_vec);

      double uni_dist;
      Utils::Vector3d uni_vec;
      uni.calculate_dist(pos, uni_dist, uni_vec);

      BOOST_CHECK_CLOSE(std::min(wall1_dist, wall2_dist), uni_dist, 1e-3);
    };

    check_union({1.2, 2.3, 4.5});
    check_union({1.2, 2.3, 5.0});
    check_union({1.2, 2.3, 5.5});
  }
}
