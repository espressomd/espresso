/*
 * Copyright (C) 2018-2022 The ESPResSo project
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

/* Unit test for Utils tuple algorithms. */

#define BOOST_TEST_MODULE Utils::tuple_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/tuple.hpp"

#include <array>
#include <functional>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

BOOST_AUTO_TEST_CASE(for_each_) {
  using Utils::for_each;

  /* l-value reference tuple */
  {
    auto a = std::array<int, 3>{{2, 3, 5}};

    for_each(
        [i = 0, a](int &e) mutable {
          a[i] = e;
          e = i++;
        },
        a);

    BOOST_CHECK_EQUAL(a[0], 0);
    BOOST_CHECK_EQUAL(a[1], 1);
    BOOST_CHECK_EQUAL(a[2], 2);
  }

  /* r-value reference */
  {
    auto t = std::make_tuple(std::make_unique<int>(5));

    for_each([](auto &e) { BOOST_CHECK_EQUAL(*e, 5); }, std::move(t));
  }

  /* move-only functor */
  {
    for_each(
        [u = std::make_unique<int>()](auto &e) { BOOST_CHECK_EQUAL(3, e); },
        std::make_pair(3, 3));
  }

  /* empty */
  {
    for_each([]() { BOOST_CHECK(false); }, std::make_tuple());
  }
}
