/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE mask test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/get.hpp>
#include <utils/mask.hpp>

#include <cstdint>
#include <string>
#include <tuple>
#include <type_traits>

BOOST_AUTO_TEST_CASE(mask_) {
  const uint8_t mask = 1u | 4u;

  auto const a = std::make_tuple(std::string("a"), 3, 4.5);
  using input_type = decltype(a);

  auto const result = Utils::mask(mask, a);
  using result_type = decltype(result);

  static_assert(std::is_same_v<input_type, result_type>);

  BOOST_CHECK_EQUAL(Utils::get<0>(result), Utils::get<0>(a));
  BOOST_CHECK_EQUAL(Utils::get<1>(result), 0);
  BOOST_CHECK_EQUAL(Utils::get<2>(result), Utils::get<2>(a));
}
