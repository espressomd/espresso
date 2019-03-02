/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::sgn test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/u32_to_u64.hpp"

BOOST_AUTO_TEST_CASE(u32_to_u64) {
  constexpr const uint64_t expected = (4ul << 32) | (11ul);
  BOOST_CHECK(expected == Utils::u32_to_u64(4u, 11u));
}

BOOST_AUTO_TEST_CASE(u64_to_u32) {
  constexpr const uint64_t expected_low = 11u;
  constexpr const uint64_t expected_high = 4u;
  constexpr const uint64_t val = (expected_high << 32) | (expected_low);

  BOOST_CHECK((std::pair<uint32_t, uint32_t>{expected_high, expected_low} ==
               Utils::u64_to_u32(val)));
}
