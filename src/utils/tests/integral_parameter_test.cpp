/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE integral_parameter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/integral_parameter.hpp>

#include <utility>

template <size_t I> struct F {
  template <class T> auto operator()(T arg) const {
    return std::make_pair(I, arg);
  }
};

BOOST_AUTO_TEST_CASE(integral_parameter_) {
  static_assert(
      std::is_same<decltype(Utils::integral_parameter<F, 1, 5>(5, 13)),
                   std::pair<size_t, int>>::value,
      "");

  BOOST_CHECK(std::make_pair(1ul, 13) ==
              (Utils::integral_parameter<F, 1, 5>(1, 13)));
  BOOST_CHECK(std::make_pair(3ul, 13) ==
              (Utils::integral_parameter<F, 1, 5>(3, 13)));
  BOOST_CHECK(std::make_pair(5ul, 13) ==
              (Utils::integral_parameter<F, 1, 5>(5, 13)));
  BOOST_CHECK_THROW((Utils::integral_parameter<F, 1, 5>(6, 13)),
                    std::exception);
}