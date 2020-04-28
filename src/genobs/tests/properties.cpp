/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <genobs/observable.hpp>
#include <genobs/properties.hpp>

#include "mock.hpp"

/* Check that the property functors correctly map to
 * the default traits */
BOOST_AUTO_TEST_CASE(properties_) {
  using namespace GenObs;

  Testing::Particle p;

  BOOST_CHECK_EQUAL(Position{}(p), traits<Testing::Particle>{}.position(p));
  BOOST_CHECK_EQUAL(Velocity{}(p), traits<Testing::Particle>{}.velocity(p));
  BOOST_CHECK_EQUAL(Mass{}(p), traits<Testing::Particle>{}.mass(p));
  BOOST_CHECK_EQUAL(Charge{}(p), traits<Testing::Particle>{}.charge(p));
  BOOST_CHECK_EQUAL(Force{}(p), traits<Testing::Particle>{}.force(p));
  BOOST_CHECK_EQUAL(DipoleMoment{}(p),
                    traits<Testing::Particle>{}.dipole_moment(p));
}
