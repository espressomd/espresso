/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

/* Unit tests for the Utils::Vector class. */

#define BOOST_TEST_MODULE Vector test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>

#include <boost/range/numeric.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

using Utils::Vector;

/* Number of nontrivial Baxter permutations of length 2n-1. (A001185) */
#define TEST_NUMBERS                                                           \
  { 0, 1, 1, 7, 21, 112, 456, 2603, 13203 }

constexpr int test_numbers[] = TEST_NUMBERS;
constexpr std::size_t n_test_numbers = sizeof(test_numbers) / sizeof(int);

template <std::size_t n> bool norm2() {
  Vector<int, n> v(std::begin(test_numbers), test_numbers + n);

  return v.norm2() == std::inner_product(v.begin(), v.end(), v.begin(), 0);
}

BOOST_AUTO_TEST_CASE(initializer_list_constructor) {
  Vector<int, n_test_numbers> v(TEST_NUMBERS);

  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));

  BOOST_CHECK_THROW(Utils::Vector2d({1., 2., 3.}), std::length_error);
}

BOOST_AUTO_TEST_CASE(iterator_constructor) {
  Vector<int, n_test_numbers> v(std::begin(test_numbers),
                                std::end(test_numbers));
  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));

  BOOST_CHECK_THROW(
      Utils::Vector2d(std::begin(test_numbers), std::end(test_numbers)),
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

BOOST_AUTO_TEST_CASE(range_constructor_test) {
  std::vector<int> test_values{TEST_NUMBERS};

  const Vector<int, n_test_numbers> v(test_values);
  BOOST_CHECK(std::equal(v.begin(), v.end(), test_numbers));
}

BOOST_AUTO_TEST_CASE(unit_vector_test) {
  BOOST_CHECK((Utils::unit_vector<int>(2) == Utils::Vector3i{0, 0, 1}));
  BOOST_CHECK_THROW(Utils::unit_vector<double>(3), std::domain_error);
  BOOST_CHECK_THROW(Utils::unit_vector<float>(-1), std::domain_error);
}

BOOST_AUTO_TEST_CASE(test_norm2) {
  BOOST_CHECK(norm2<1>());
  BOOST_CHECK(norm2<2>());
  BOOST_CHECK(norm2<3>());
  BOOST_CHECK(norm2<4>());
}

BOOST_AUTO_TEST_CASE(normalize) {
  Utils::Vector3d v{1, 2, 3};
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
  Utils::Vector3i v1{1, 2, 3};
  Utils::Vector3i v2{4, 5, 6};

  BOOST_CHECK(((v1 + v2) == Utils::Vector3i{5, 7, 9}));
  BOOST_CHECK(((v1 - v2) == Utils::Vector3i{-3, -3, -3}));
  BOOST_CHECK(((-v1) == Utils::Vector3i{-1, -2, -3}));

  /* Mixed types */
  {
    BOOST_CHECK((Utils::Vector3d{1., 2., 3.} * 4) ==
                (Utils::Vector3d{1. * 4, 2. * 4, 3. * 4}));
    BOOST_CHECK((Utils::Vector3d{1., 2., 3.} + Utils::Vector3i{11, 12, 13}) ==
                (Utils::Vector3d{1. + 11, 2. + 12, 3. + 13}));
    BOOST_CHECK((Utils::Vector3d{1., 2., 3.} - Utils::Vector3i{11, 12, 13}) ==
                (Utils::Vector3d{1. - 11, 2. - 12, 3. - 13}));
  }

  {
    auto v3 = v1;
    BOOST_CHECK((v1 + v2) == (v3 += v2));
  }

  {
    auto v3 = v1;
    BOOST_CHECK((v1 - v2) == (v3 -= v2));
  }

  BOOST_CHECK(((2 * v1) == Utils::Vector3i{2, 4, 6}));
  BOOST_CHECK(((v1 * 2) == Utils::Vector3i{2, 4, 6}));

  {
    Utils::Vector3i v3{2, 4, 6};
    auto v4 = 2 * v3;
    BOOST_CHECK(v4 == (v3 *= 2));
  }

  {
    Utils::Vector3i v3{2, 4, 6};
    auto v4 = v3 / 2;
    BOOST_CHECK(v4 == (v3 /= 2));
  }

  BOOST_CHECK((sqrt(Utils::Vector3d{1., 2., 3.}) ==
               Utils::Vector3d{sqrt(1.), sqrt(2.), sqrt(3.)}));

  /* modulo */
  {
    Utils::Vector3i v3{2, 7, 8};
    Utils::Vector3i v4{1, 2, 3};

    auto const res = v3 % v4;

    BOOST_CHECK_EQUAL(res[0], v3[0] % v4[0]);
    BOOST_CHECK_EQUAL(res[1], v3[1] % v4[1]);
    BOOST_CHECK_EQUAL(res[2], v3[2] % v4[2]);
  }
}

BOOST_AUTO_TEST_CASE(broadcast) {
  const auto fives = Vector<int, 11>::broadcast(5);

  BOOST_CHECK(fives.size() == 11);
  for (auto const &e : fives) {
    BOOST_CHECK(5 == e);
  }
}

BOOST_AUTO_TEST_CASE(swap) {
  const auto cv1 = Utils::Vector3i{1, 2, 3};
  const auto cv2 = Utils::Vector3i{4, 5, 6};

  auto v1 = cv1;
  auto v2 = cv2;

  v1.swap(v2);

  BOOST_CHECK(v1 == cv2);
  BOOST_CHECK(v2 == cv1);
}

BOOST_AUTO_TEST_CASE(decay_to_scalar_test) {
  {
    using original_t = Vector<int, 1>;
    using decayed_t = typename Utils::decay_to_scalar<original_t>::type;

    static_assert(std::is_same_v<int, decayed_t>);
  }
  BOOST_TEST_PASSPOINT();
  {
    using original_t = Utils::Vector3i;
    using decayed_t = typename Utils::decay_to_scalar<original_t>::type;

    static_assert(std::is_same_v<original_t, decayed_t>);
  }
}

BOOST_AUTO_TEST_CASE(vector_broadcast) {
  Utils::Vector2d const v = Utils::Vector2d::broadcast(1.4);

  BOOST_CHECK_EQUAL(v[0], 1.4);
  BOOST_CHECK_EQUAL(v[1], 1.4);
}

BOOST_AUTO_TEST_CASE(products) {
  /* Types */
  {
    static_assert(
        std::is_same_v<decltype(Utils::Vector3d{} * Utils::Vector3d{}),
                       double>);
    static_assert(
        std::is_same_v<decltype(Utils::Vector3d{} * Utils::Vector3i{}),
                       double>);
    static_assert(
        std::is_same_v<decltype(Vector<std::complex<float>, 2>{} * 3.f),
                       Vector<std::complex<float>, 2>>);
  }

  /* Vector-Vector */
  {
    auto const v1 = Utils::Vector3d{1., 2., 3.};
    auto const v2 = Utils::Vector3d{4.1, 5.2, 6.3};
    auto const v3 = Utils::Vector3i{11, 12, 13};
    BOOST_CHECK_EQUAL(v1 * v2, boost::inner_product(v1, v2, 0.));
    BOOST_CHECK_EQUAL(v1 * v3, boost::inner_product(v1, v3, 0.));
  }
}

BOOST_AUTO_TEST_CASE(conversion) {
  using Utils::Vector3d;
  using Utils::Vector3f;

  /* The Vector class has a user-defined conversion method,
   * however static_cast() will not use it. Instead, the
   * range-based constructor is called. */

  auto const orig = Vector3d{1., 2., 3.};
  auto const expected =
      Vector3f{static_cast<float>(orig[0]), static_cast<float>(orig[1]),
               static_cast<float>(orig[2])};

  // check range-based conversion
  {
    auto const result = static_cast<Vector3f>(orig);
    BOOST_TEST(result == expected);
  }

  // check cast operator
  {
    auto const result = orig.operator Utils::Vector3f();
    BOOST_TEST(result == expected);
  }

  // check vector conversion
  {
    auto const result = Utils::Vector3d{static_cast<std::vector<double>>(orig)};
    BOOST_TEST(result == orig);
  }

  // check vector conversion
  {
    auto const result = Utils::Vector3d{orig.as_vector()};
    BOOST_TEST(result == orig);
  }
}

BOOST_AUTO_TEST_CASE(vector_product_test) {
  auto const v1 = Utils::Vector3d{1., 2., 3.};
  auto const v2 = Utils::Vector3d{4.1, 5.2, 6.3};

  auto const res = vector_product(v1, v2);
  BOOST_CHECK_SMALL(std::abs(v1 * res), 1e-14);
  BOOST_CHECK_SMALL(std::abs(v2 * res), 1e-14);
}

BOOST_AUTO_TEST_CASE(product_test) {
  auto const v1 = Utils::Vector<int, 2>{8, -9};
  auto const v2 = Utils::Vector<int, 2>{0, 6};

  BOOST_CHECK_EQUAL(Utils::product(v1), -8 * 9);
  BOOST_CHECK_EQUAL(Utils::product(v2), 0);
}

BOOST_AUTO_TEST_CASE(hadamard_product_test) {
  auto const v1 = Utils::Vector<int, 2>{8, 9};
  auto const v2 = Utils::Vector<int, 2>{5, 6};
  auto const a = 6;

  auto const res1 = Utils::hadamard_product(v1, v2);
  BOOST_CHECK_EQUAL(res1[0], v1[0] * v2[0]);
  BOOST_CHECK_EQUAL(res1[1], v1[1] * v2[1]);

  auto const res2 = Utils::hadamard_product(a, v1);
  BOOST_CHECK_EQUAL(res2[0], a * v1[0]);
  BOOST_CHECK_EQUAL(res2[1], a * v1[1]);

  auto const res3 = Utils::hadamard_product(v1, a);
  BOOST_CHECK_EQUAL(res3[0], a * v1[0]);
  BOOST_CHECK_EQUAL(res3[1], a * v1[1]);

  auto const res4 = Utils::hadamard_product(a, a);
  BOOST_CHECK_EQUAL(res4, a * a);
}

BOOST_AUTO_TEST_CASE(hadamard_division_test) {
  auto const v1 = Utils::Vector<int, 2>{16, 32};
  auto const v2 = Utils::Vector<int, 2>{8, 4};
  auto const a = 2;

  auto const res1 = Utils::hadamard_division(v1, v2);
  BOOST_CHECK_EQUAL(res1[0], v1[0] / v2[0]);
  BOOST_CHECK_EQUAL(res1[1], v1[1] / v2[1]);

  auto const res2 = Utils::hadamard_division(a, v1);
  BOOST_CHECK_EQUAL(res2[0], a / v1[0]);
  BOOST_CHECK_EQUAL(res2[1], a / v1[1]);

  auto const res3 = Utils::hadamard_division(v1, a);
  BOOST_CHECK_EQUAL(res3[0], v1[0] / a);
  BOOST_CHECK_EQUAL(res3[1], v1[1] / a);

  auto const res4 = Utils::hadamard_division(a, a);
  BOOST_CHECK_EQUAL(res4, 1);
}

BOOST_AUTO_TEST_CASE(type_deduction) {
  static_assert(
      std::is_same_v<
          typename boost::qvm::deduce_vec<Utils::Vector<double, 3>, 3>::type,
          Utils::Vector<double, 3>>);
  BOOST_TEST_PASSPOINT();
}
