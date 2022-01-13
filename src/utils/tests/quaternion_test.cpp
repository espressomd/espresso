/*
 * Copyright (C) 2019 The ESPResSo project
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

#define BOOST_TEST_MODULE quaternion test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/quaternion.hpp"
#include "utils/quaternion.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <limits>

constexpr int w = -5;
constexpr int x = 2;
constexpr int y = 3;
constexpr int z = 4;

Utils::Quaternion<int> scalar_quat{w, 0, 0, 0};
Utils::Quaternion<int> full_quat{w, x, y, z};
Utils::Quaternion<int> vector_quat{0, x, y, z};

BOOST_AUTO_TEST_CASE(multiply_quaternions) {
  /* identities */
  BOOST_CHECK(Utils::Quaternion<int>::identity() * scalar_quat == scalar_quat);
  BOOST_CHECK(Utils::Quaternion<int>::identity() * full_quat == full_quat);
  BOOST_CHECK(scalar_quat * scalar_quat == w * scalar_quat);
  BOOST_CHECK(scalar_quat * vector_quat == w * vector_quat);
  BOOST_CHECK(Utils::Quaternion<int>::zero() * full_quat ==
              Utils::Quaternion<int>::zero());
  BOOST_CHECK(vector_quat * vector_quat ==
              -vector_quat.norm2() * Utils::Quaternion<int>::identity());

  /* other */
  Utils::Quaternion<int> const reference_quat{{-4, -20, -30, -40}};
  BOOST_CHECK(full_quat * full_quat == reference_quat);
}

BOOST_AUTO_TEST_CASE(convert_quaternion_to_director) {
  using Utils::convert_quaternion_to_director;
  /* identities */
  BOOST_CHECK(
      convert_quaternion_to_director(Utils::Quaternion<int>::identity()) ==
      (Utils::Vector3i{{0, 0, 1}}));
  BOOST_CHECK(convert_quaternion_to_director(scalar_quat) ==
              (Utils::Vector3i{{0, 0, w * w}}));
  BOOST_CHECK(
      convert_quaternion_to_director(vector_quat) ==
      (Utils::Vector3i{{2 * x * z, 2 * y * z, -x * x - y * y + z * z}}));
  /* other */
  Utils::Vector3i const reference_director{{-14, 44, 28}};
  BOOST_CHECK(convert_quaternion_to_director(full_quat) == reference_director);
}

BOOST_AUTO_TEST_CASE(convert_director_to_quaternion) {
  using Quat = Utils::Quaternion<double>;
  using Utils::convert_director_to_quaternion;
  using Utils::Vector3d;
  double const cos_pi_4 = std::sqrt(2.) / 2.;
  constexpr double eps = std::numeric_limits<double>::epsilon();
#define CHECK_QUAT(input, ref)                                                 \
  BOOST_CHECK_LE((convert_director_to_quaternion(input) - (ref)).norm2(), eps);
  /* identities */
  CHECK_QUAT((Vector3d{{0, 0, 0}}), (Quat{{1, 0, 0, 0}}));
  CHECK_QUAT((Vector3d{{0, 0, +1}}), (Quat{{1, 0, 0, 0}}));
  CHECK_QUAT((Vector3d{{0, 0, -1}}), (Quat{{0, -1, 0, 0}}));
  CHECK_QUAT((Vector3d{{+1, 0, 0}}), (Quat{{+1, -1, +1, -1}} / 2.));
  CHECK_QUAT((Vector3d{{-1, 0, 0}}), (Quat{{-1, +1, +1, -1}} / 2.));
  CHECK_QUAT((Vector3d{{0, +1, 0}}), (Quat{{+1, -1, 0, 0}} * cos_pi_4));
  CHECK_QUAT((Vector3d{{0, -1, 0}}), (Quat{{0, 0, +1, -1}} * cos_pi_4));
  /* self-consistency */
  using Utils::convert_quaternion_to_director;
  for (int i = -2; i < 3; ++i) {
    for (int j = -2; j < 3; ++j) {
      for (int k = -2; k < 3; ++k) {
        if (i == 0 and j == 0 and k == 0) {
          continue;
        }
        // here the case where j is zero is also tested
        Vector3d input{{static_cast<double>(i), static_cast<double>(j),
                        static_cast<double>(k)}};
        Vector3d output = convert_quaternion_to_director(
            convert_director_to_quaternion(input));
        input.normalize();
        output.normalize();
        BOOST_CHECK_LE((input - output).norm2(), eps);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(quat_type) {
  Utils::Quaternion<double> test{1, 2, 3, 4};
  BOOST_CHECK(test[0] == 1);
  test.normalize();
  BOOST_CHECK_LE(test.norm() - 1.0, std::numeric_limits<double>::epsilon());
  BOOST_CHECK((Utils::Quaternion<int>::identity() ==
               Utils::Quaternion<int>{1, 0, 0, 0}));
  BOOST_CHECK(
      (Utils::Quaternion<int>::zero() == Utils::Quaternion<int>{0, 0, 0, 0}));
  BOOST_CHECK((Utils::Quaternion<double>{1, 0, 0, 0} ==
               Utils::Quaternion<double>{2, 0, 0, 0}.normalized()));
  BOOST_CHECK_SMALL(
      (Utils::Quaternion<double>{2, 1, 3, 4}.normalized().norm() - 1.0),
      std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(type_deduction) {
  using Utils::Quaternion;
  using Utils::Vector;
  static_assert(
      std::is_same<typename boost::qvm::deduce_scalar<Quaternion<double>,
                                                      Quaternion<double>>::type,
                   double>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_scalar<Quaternion<double>,
                                                      Quaternion<int>>::type,
                   double>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_scalar<Quaternion<double>,
                                                      Vector<double, 3>>::type,
                   double>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_vec2<Quaternion<double>,
                                                    Vector<double, 3>, 3>::type,
                   Vector<double, 3>>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_vec2<Quaternion<double>,
                                                    Vector<int, 3>, 3>::type,
                   Vector<double, 3>>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_quat<Quaternion<double>>::type,
                   Quaternion<double>>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_quat2<Quaternion<double>,
                                                     Quaternion<double>>::type,
                   Quaternion<double>>::value,
      "");
  static_assert(
      std::is_same<typename boost::qvm::deduce_quat2<Quaternion<double>,
                                                     Quaternion<int>>::type,
                   Quaternion<double>>::value,
      "");
  BOOST_TEST_PASSPOINT();
}
