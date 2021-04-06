/*
 * Copyright (C) 2021 The ESPResSo project
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

#define BOOST_TEST_MODULE CylindricalTransformationParameters test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/math/cylindrical_transformation_parameters.hpp>

#include <stdexcept>

BOOST_AUTO_TEST_CASE(getters_test) {
  using CTP = Utils::CylindricalTransformationParameters;
  Utils::Vector3d const center{-1., 0., 1.};
  Utils::Vector3d const axis{1., 0., 0.};
  Utils::Vector3d const orientation{0., -1., 0.};
  auto const ctp = CTP(center, axis, orientation);
  BOOST_TEST(ctp.center() == center, boost::test_tools::per_element());
  BOOST_TEST(ctp.axis() == axis, boost::test_tools::per_element());
  BOOST_TEST(ctp.orientation() == orientation,
             boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(exceptions_test) {
  using CTP = Utils::CylindricalTransformationParameters;
  BOOST_CHECK_THROW(CTP({0, 0, 0}, {1, 0, 0}, {1, 0, 0}), std::runtime_error);
  BOOST_CHECK_THROW(CTP({0, 0, 0}, {2, 0, 0}, {0, 1, 0}), std::runtime_error);
  BOOST_CHECK_THROW(CTP({0, 0, 0}, {1, 0, 0}, {0, 2, 0}), std::runtime_error);
}
