/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::AS_erfc_part test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/AS_erfc_part.hpp"
using Utils::AS_erfc_part;

#include <cmath>

/* Check that it can be used in constexpr context */
static_assert((AS_erfc_part(0.), true), "");

BOOST_AUTO_TEST_CASE(approx) {
  for (double x = 0.0; x <= 1.; x += 0.01) {
    auto const approx = AS_erfc_part(x);
    auto const exact = std::exp(x * x) * std::erfc(x);
    BOOST_CHECK(std::abs(approx - exact) < 5.e-7);
  }
}
