/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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

/** \file
 * Unit tests for the Utils::Vector class.
 *
 */

#define BOOST_TEST_MODULE Vector test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

#include "utils/Vector.hpp"

/** Number of nontrivial Baxter permutations of length 2n-1. (A001185) */
#define TEST_NUMBERS                                                           \
  { 0, 1, 1, 7, 21, 112, 456, 2603, 13203 }
#define TEST_NUMBERS_PARTIAL_NORM2                                             \
  { 0, 1, 2, 51, 492, 13036 }
constexpr int test_numbers[] = TEST_NUMBERS;
constexpr int n_test_numbers = sizeof(test_numbers) / sizeof(int);

template <int n> bool norm2() {
  Vector<int, n> v(std::begin(test_numbers), test_numbers + n);

  return v.norm2() == std::inner_product(v.begin(), v.end(), v.begin(), 0);
}

BOOST_AUTO_TEST_CASE(initializer_list_constructor) {
  Vector<int, n_test_numbers> v(TEST_NUMBERS);

  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));

  BOOST_CHECK_THROW(Vector2d({1., 2., 3.}), std::length_error);
}

BOOST_AUTO_TEST_CASE(iterator_constructor) {
  Vector<int, n_test_numbers> v(std::begin(test_numbers),
                                std::end(test_numbers));
  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));

  BOOST_CHECK_THROW(Vector2d(std::begin(test_numbers), std::end(test_numbers)),
                    std::length_error);
}

BOOST_AUTO_TEST_CASE(const_iterator_constructor) {
  /* {begin,end}() const variant */
  {
    const Vector<int, n_test_numbers> v(std::begin(test_numbers),
                                        std::end(test_numbers));
    BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));
  }

  /* {cbegin,cend}() const variant */
  {
    Vector<int, n_test_numbers> v(std::begin(test_numbers),
                                  std::end(test_numbers));
    BOOST_CHECK(std::equal(v.cbegin(), v.cend(), test_numbers));
  }
}

BOOST_AUTO_TEST_CASE(default_constructor_test) {
  Vector<int, 1> v2;
  BOOST_CHECK(v2.size() == 1);
  Vector<int, 2> v3;
  BOOST_CHECK(v3.size() == 2);
  Vector<int, 11> v4;
  BOOST_CHECK(v4.size() == 11);
}

BOOST_AUTO_TEST_CASE(test_norm2) {
  BOOST_CHECK(norm2<1>());
  BOOST_CHECK(norm2<2>());
  BOOST_CHECK(norm2<3>());
  BOOST_CHECK(norm2<4>());
}

BOOST_AUTO_TEST_CASE(normalize) {
  Vector3d v{1, 2, 3};
  v.normalize();

  BOOST_CHECK((v.norm2() - 1.0) <= std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(comparison_operators) {
  Vector<int, 5> v1{1, 2, 3, 4, 5};
  Vector<int, 5> v2{6, 7, 8, 9, 10};
  Vector<int, 5> v3{1, 2, 3, 5, 6};

  BOOST_CHECK(v1 < v2);
  BOOST_CHECK(!(v1 < v1));
  BOOST_CHECK(v1 <= v2);
  BOOST_CHECK(v1 <= v1);
  BOOST_CHECK(v2 > v1);
  BOOST_CHECK(!(v2 > v2));
  BOOST_CHECK(v2 >= v1);
  BOOST_CHECK(v2 >= v2);
  BOOST_CHECK(v1 != v2);
  BOOST_CHECK(v1 != v3);
  BOOST_CHECK(v1 <= v3);
  BOOST_CHECK(!(v1 == v2));
  BOOST_CHECK(v1 == v1);
}

BOOST_AUTO_TEST_CASE(algebraic_operators) {
  Vector3i v1{1, 2, 3};
  Vector3i v2{4, 5, 6};

  BOOST_CHECK((v1 * v2) ==
              std::inner_product(v1.begin(), v1.end(), v2.begin(), 0));

  BOOST_CHECK(((v1 + v2) == Vector3i{5, 7, 9}));
  BOOST_CHECK(((v1 - v2) == Vector3i{-3, -3, -3}));
  BOOST_CHECK(((-v1) == Vector3i{-1, -2, -3}));

  {
    auto v3 = v1;
    BOOST_CHECK((v1 + v2) == (v3 += v2));
  }

  {
    auto v3 = v1;
    BOOST_CHECK((v1 - v2) == (v3 -= v2));
  }

  BOOST_CHECK(((2 * v1) == Vector3i{2, 4, 6}));
  BOOST_CHECK(((v1 * 2) == Vector3i{2, 4, 6}));

  {
    Vector3i v1{2, 4, 6};
    auto v2 = 2 * v1;
    BOOST_CHECK(v2 == (v1 *= 2));
  }

  {
    Vector3i v1{2, 4, 6};
    auto v2 = v1 / 2;
    BOOST_CHECK(v2 == (v1 /= 2));
  }

  BOOST_CHECK(
      (sqrt(Vector3d{1., 2., 3.}) == Vector3d{sqrt(1.), sqrt(2.), sqrt(3.)}));
}

BOOST_AUTO_TEST_CASE(broadcast) {
  const auto fives = Vector<int, 11>::broadcast(5);

  BOOST_CHECK(fives.size() == 11);
  for (auto const &e : fives) {
    BOOST_CHECK(5 == e);
  }
}

BOOST_AUTO_TEST_CASE(swap) {
  const auto cv1 = Vector3i{1, 2, 3};
  const auto cv2 = Vector3i{4, 5, 6};

  auto v1 = cv1;
  auto v2 = cv2;

  v1.swap(v2);

  BOOST_CHECK(v1 == cv2);
  BOOST_CHECK(v2 == cv1);
}

BOOST_AUTO_TEST_CASE(decay_to_scalar_test) {
  {
    using original_t = Vector<int, 1>;
    using decayed_t = typename decay_to_scalar<original_t>::type;

    static_assert(std::is_same<int, decayed_t>::value, "");
  }

  {
    using original_t = Vector3i;
    using decayed_t = typename decay_to_scalar<original_t>::type;

    static_assert(std::is_same<original_t, decayed_t>::value, "");
  }
}

BOOST_AUTO_TEST_CASE(vector_broadcast) {
  Vector2d const v = Vector2d::broadcast(1.4);

  BOOST_CHECK_EQUAL(v[0], 1.4);
  BOOST_CHECK_EQUAL(v[1], 1.4);
}
