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

#define BOOST_TEST_MODULE Utils::Array test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Array.hpp>
#include <utils/get.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>

using Utils::Array;

BOOST_AUTO_TEST_CASE(const_expr_ctor) {
  static_assert(4 == Array<int, 4>().size());
  static_assert(4 == Array<int, 4>().max_size());
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(array_ctor) {
  Array<int, 4> a;
  Array<int, 0> b;

  BOOST_CHECK_EQUAL(a.size(), 4);
  BOOST_CHECK_EQUAL(a.max_size(), 4);
  BOOST_CHECK_EQUAL(b.size(), 0);
  BOOST_CHECK_EQUAL(b.max_size(), 0);
}

BOOST_AUTO_TEST_CASE(iterators) {
  auto a = Array<int, 4>{{{1, 2, 3, 4}}};

  BOOST_CHECK_EQUAL(*(a.begin()), 1);
  BOOST_CHECK_EQUAL(*(a.cbegin()), 1);
  BOOST_CHECK_EQUAL(*(a.end() - 1), 4);
  BOOST_CHECK_EQUAL(*(a.cend() - 1), 4);
}

BOOST_AUTO_TEST_CASE(element_access) {
  auto a = Array<int, 5>{{{5, 6, 7, 8, 9}}};
  auto const &b = a;

  int c = 5;
  int j = 0;
  for (int i : a) {
    BOOST_CHECK_EQUAL(i, c);
    BOOST_CHECK_EQUAL(a[j], c);
    BOOST_CHECK_EQUAL(b[j], c);
    BOOST_CHECK_EQUAL(a.at(j), c);
    BOOST_CHECK_EQUAL(b.at(j), c);
    ++c;
    ++j;
  }

  BOOST_CHECK_THROW(a.at(a.size()), std::out_of_range);
  BOOST_CHECK_THROW(b.at(b.size()), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(fill) {
  Array<int, 10> a{};
  a.fill(10);
  for (int i : a) {
    BOOST_CHECK_EQUAL(i, 10);
  }
}

BOOST_AUTO_TEST_CASE(broadcast) {
  constexpr auto a = Array<int, 3>::broadcast(5);
  static_assert(a[0] == 5);
  static_assert(a[1] == 5);
  static_assert(a[2] == 5);
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(serialization) {
  Array<double, 24> a;
  std::iota(std::begin(a), std::end(a), 0);

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << a;

  boost::archive::text_iarchive in_ar(stream);
  decltype(a) b;
  in_ar >> b;

  BOOST_CHECK(std::equal(std::begin(a), std::end(a), std::begin(b)));
}

BOOST_AUTO_TEST_CASE(zero_size) {
  Array<int, 0> const a{};
  BOOST_CHECK(a.empty());
}

BOOST_AUTO_TEST_CASE(tuple_protocol) {
  using A = Utils::Array<int, 4>;

  static_assert(std::is_same_v<Utils::tuple_element_t<0, A>, int>);
  static_assert(std::is_same_v<Utils::tuple_element_t<1, A>, int>);
  static_assert(A{}.size() == Utils::tuple_size<A>::value);

  BOOST_CHECK_EQUAL(Utils::get<1>(A{{{1, 2, 3, 4}}}), 2);
}

BOOST_AUTO_TEST_CASE(streaming_operator) {
  {
    auto const a = Utils::Array<int, 1>{{{1}}};

    std::stringstream ss;
    ss << a;

    BOOST_CHECK_EQUAL(ss.str(), "1");
  }

  {
    auto const a = Utils::Array<int, 3>{{{1, 2, 3}}};

    std::stringstream ss;
    ss << a;

    BOOST_CHECK_EQUAL(ss.str(), "1, 2, 3");
  }
}

BOOST_AUTO_TEST_CASE(formatter_and_streaming_operator) {
  {
    auto const a = Utils::Array<int, 1>{{{1}}};

    std::stringstream ss;
    ss << a.formatter("xyz") << a;

    BOOST_CHECK_EQUAL(ss.str(), "1");
  }

  {
    auto const a = Utils::Array<int, 3>{{{1, 2, 3}}};

    std::stringstream ss;
    ss << a.formatter(" + ") << a;

    BOOST_CHECK_EQUAL(ss.str(), "1 + 2 + 3");
  }
}
