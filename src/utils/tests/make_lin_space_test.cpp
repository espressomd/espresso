/*
Copyright (C) 2019 The ESPResSo project

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

#define BOOST_TEST_MODULE make_lin_space test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <cmath>

#include <utils/math/make_lin_space.hpp>

#include <vector>

BOOST_AUTO_TEST_CASE(make_lin_space_test) {
  using Utils::make_lin_space;

  /* With endpoint */
  {
    auto const start = 1.;
    auto const stop = 2.;
    auto const num = 13;

    auto const lin_space =
        make_lin_space(start, stop, num, /* endpoint */ true);
    BOOST_CHECK_EQUAL(lin_space.size(), num);

    std::vector<double> values(lin_space.begin(), lin_space.end());
    BOOST_CHECK_EQUAL(values.front(), start);
    BOOST_CHECK_EQUAL(values.back(), stop);

    auto const dx = (stop - start) / (num - 1);
    for (int i = 0; i < values.size(); i++) {
      BOOST_CHECK(std::fabs(start + i * dx - values.at(i)) <=
                  std::numeric_limits<double>::epsilon());
    }
  }

  /* Without endpoint */
  {
    auto const start = 1.;
    auto const stop = 2.;
    auto const num = 13;

    auto const lin_space =
        make_lin_space(start, stop, num, /* endpoint */ false);
    BOOST_CHECK_EQUAL(lin_space.size(), num);

    std::vector<double> values(lin_space.begin(), lin_space.end());
    BOOST_CHECK_EQUAL(values.front(), start);
    BOOST_CHECK_LT(values.back(), stop);

    auto const dx = (stop - start) / num;
    for (int i = 0; i < values.size(); i++) {
      BOOST_CHECK(std::fabs(start + i * dx - values.at(i)) <=
                  std::numeric_limits<double>::epsilon());
    }
  }
}
