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

#define BOOST_TEST_MODULE Utils::Span test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Span.hpp"
using Utils::Span;

#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

BOOST_AUTO_TEST_CASE(const_expr_ctor) {
  static_assert(4 == Span<int>(nullptr, 4).size());
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(array_ctor) {
  BOOST_CHECK((std::is_constructible_v<Span<const int>, int[3]>));
  BOOST_CHECK((std::is_constructible_v<Span<const int>, const int[3]>));
  BOOST_CHECK(not(std::is_constructible_v<Span<int>, const int[3]>));
  BOOST_CHECK((std::is_convertible_v<int[3], Span<const int>>));
  BOOST_CHECK((std::is_convertible_v<const int[3], Span<const int>>));

  int a[4] = {1, 2, 3, 4};
  Span<int> s(a);

  BOOST_CHECK_EQUAL(s.data(), a);
  BOOST_CHECK_EQUAL(s.size(), 4);
}

BOOST_AUTO_TEST_CASE(ctor) {
  /* Container conversion rules */
  {
    BOOST_CHECK((std::is_constructible_v<Span<const int>, std::vector<int>>));
    BOOST_CHECK(
        (std::is_constructible_v<Span<const int>, const std::vector<int>>));
    BOOST_CHECK(
        not(std::is_constructible_v<Span<int>, const std::vector<int>>));
    BOOST_CHECK((std::is_convertible_v<std::vector<int>, Span<const int>>));
    BOOST_CHECK(
        (std::is_convertible_v<const std::vector<int>, Span<const int>>));
  }

  /* from ptr + size */
  {
    std::vector<int> v(23);

    auto s = Span<int>(v.data(), v.size());

    BOOST_CHECK(v.size() == s.size());
    BOOST_CHECK(v.data() == s.data());
  }

  /* From container */
  {
    std::vector<int> v{{1, 2, 3}};
    auto s = Span<int>(v);

    BOOST_CHECK(v.size() == s.size());
    BOOST_CHECK(v.data() == s.data());
  }
}

BOOST_AUTO_TEST_CASE(iterators) {
  int dummy;
  auto const p = &dummy;
  auto const size = 11u;

  auto s = Span<int>(p, size);

  BOOST_CHECK(s.begin() == p);
  BOOST_CHECK(s.cbegin() == p);
  BOOST_CHECK(s.end() == (p + size));
  BOOST_CHECK(s.cend() == (p + size));
}

BOOST_AUTO_TEST_CASE(element_access) {
  std::vector<int> v(11);
  std::iota(v.begin(), v.end(), 5);

  auto s = Span<int>(v.data(), v.size());

  for (Span<int>::size_type i = 0; i < s.size(); ++i) {
    BOOST_CHECK(v.at(i) == s[i]);
    BOOST_CHECK(v.at(i) == s.at(i));
  }

  BOOST_CHECK_THROW(s.at(s.size()), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(make_span_) {
  using std::declval;
  using Utils::make_span;

  static_assert(std::is_same_v<decltype(make_span(declval<int *>(),
                                                  declval<std::size_t>())),
                               Span<int>>);
  static_assert(std::is_same_v<decltype(make_span(declval<const int *>(),
                                                  declval<std::size_t>())),
                               Span<const int>>);

  /* From pointer and size */
  {
    const int p = 5;
    auto s = make_span(&p, 1);
    BOOST_CHECK_EQUAL(s.data(), &p);
    BOOST_CHECK_EQUAL(s.size(), 1);
  }

  /* From container */
  {
    std::vector<int> vec(5);
    auto result = make_span(vec);
    BOOST_CHECK_EQUAL(result.data(), vec.data());
    BOOST_CHECK_EQUAL(result.size(), vec.size());
  }
}

BOOST_AUTO_TEST_CASE(make_const_span_) {
  using std::declval;
  using Utils::make_const_span;

  static_assert(std::is_same_v<decltype(make_const_span(
                                   declval<int *>(), declval<std::size_t>())),
                               Span<const int>>);
  static_assert(
      std::is_same_v<decltype(make_const_span(declval<const int *>(),
                                              declval<std::size_t>())),
                     Span<const int>>);

  {
    const int p = 5;
    auto s = make_const_span(&p, 1);
    BOOST_CHECK_EQUAL(s.data(), &p);
    BOOST_CHECK_EQUAL(s.size(), 1);
  }

  {
    std::vector<int> vec(5);
    auto result = make_const_span(vec);
    BOOST_CHECK_EQUAL(result.data(), vec.data());
    BOOST_CHECK_EQUAL(result.size(), vec.size());
  }
}
