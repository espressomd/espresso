/*
 * Copyright (C) 2018-2019 The ESPResSo project
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

BOOST_AUTO_TEST_CASE(for_each) {
  using Utils::for_each;

  /* l-value reference tuple */
  {
    auto a = std::array<int, 3>{2, 3, 5};

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

BOOST_AUTO_TEST_CASE(apply) {
  using Utils::apply;

  /* constexpr */
  { static_assert(apply(std::plus<>(), std::array<int, 2>{3, 8}) == 11, ""); }

  /* l-value reference */
  {
    auto t = std::make_tuple(4, 0, 7);

    apply(
        [](int &a, int &b, int &c) {
          BOOST_CHECK_EQUAL(a, 4);
          BOOST_CHECK_EQUAL(b, 0);
          BOOST_CHECK_EQUAL(c, 7);

          a = -1;
          b = -2;
          c = -3;
        },
        t);

    BOOST_CHECK_EQUAL(std::get<0>(t), -1);
    BOOST_CHECK_EQUAL(std::get<1>(t), -2);
    BOOST_CHECK_EQUAL(std::get<2>(t), -3);
  }

  /* r-value reference */
  {
    apply([](auto &&a) { BOOST_CHECK_EQUAL(*a, 4); },
          std::make_tuple(std::make_unique<int>(4

                                                )));
  }

  /* empty */
  {
    bool called = false;
    apply([&called]() { called = true; }, std::make_tuple());

    BOOST_CHECK(called);
  }
}

BOOST_AUTO_TEST_CASE(find_if_) {
  {
    auto const result = Utils::find_if([](int e) { return e == 2; },
                                       std::array<int, 4>{1, 2, 3, 4},
                                       [](int e) { BOOST_CHECK_EQUAL(e, 2); });
    BOOST_CHECK(result);
  }

  {
    auto const result = Utils::find_if([](int e) { return e == 5; },
                                       std::array<int, 4>{1, 2, 3, 4},
                                       [](int e) { BOOST_CHECK(false); });
    BOOST_CHECK(not result);
  }
}

BOOST_AUTO_TEST_CASE(filter_) {
  using Utils::filter;

  constexpr auto const expected = std::make_tuple(1, 2u);
  constexpr auto const result =
      filter<std::is_integral>(std::make_tuple(1, 1.5, 2u, 2.5));

  BOOST_CHECK(expected == result);
}