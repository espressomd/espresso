/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::Interaction::LennardJones test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <limits>

#include "utils/Vector.hpp"
#include "utils/interaction/LennardJones.hpp"

BOOST_AUTO_TEST_CASE(energy) {
  Utils::Interaction::LennardJones lennard_jones;
  lennard_jones.sigma = 1.5;
  lennard_jones.epsilon = 2.0;
  // check the minimum
  BOOST_CHECK(lennard_jones.U(std::pow(2.0, 1. / 6) * lennard_jones.sigma) ==
              -lennard_jones.epsilon);
  // check the root
  BOOST_CHECK(lennard_jones.U(lennard_jones.sigma) == 0.0);
}

BOOST_AUTO_TEST_CASE(force) {
  Utils::Interaction::LennardJones lennard_jones;
  lennard_jones.sigma = 1.5;
  lennard_jones.epsilon = 2.0;
  // check the minimum
  auto const minimum = std::pow(2.0, 1. / 6) * lennard_jones.sigma;
  auto const force =
      lennard_jones.F(minimum, Utils::Vector3d{minimum, 0.0, 0.0});
  auto const expected = Utils::Vector3d::broadcast(0.0);
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(force[i] - expected[i], 1e-10);
  }
  // Check for sign change in the force
  auto const epsilon = std::numeric_limits<double>::epsilon();
  auto const force1 = lennard_jones.F(
      minimum - epsilon, Utils::Vector3d{minimum - epsilon, 0.0, 0.0});
  auto const force2 = lennard_jones.F(
      minimum + epsilon, Utils::Vector3d{minimum + epsilon, 0.0, 0.0});
  BOOST_CHECK(force1[0] > 0.0);
  BOOST_CHECK(force2[0] < 0.0);
}
