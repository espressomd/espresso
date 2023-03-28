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

#define BOOST_TEST_MODULE Utils::Matrix test
#define BOOST_TEST_DYN_LINK

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>

BOOST_AUTO_TEST_CASE(matrix) {
  Utils::Matrix<int, 2, 2> mat2{{8, 2}, {3, 4}};
  BOOST_CHECK((mat2(0, 0) == 8));
  BOOST_CHECK((mat2(1, 0) == 3));
  BOOST_CHECK((mat2(0, 1) == 2));
  BOOST_CHECK((mat2(1, 1) == 4));

  Utils::Matrix<int, 3, 3> mat{1, 2, 3, 4, 5, 6, 7, 8, 9};

  Utils::Vector3i diagonal = boost::qvm::diag(mat);
  BOOST_CHECK((diagonal == Utils::Vector3i{1, 5, 9}));

  Utils::Vector3i row = boost::qvm::row<0>(mat);
  BOOST_CHECK((row == Utils::Vector3i{1, 2, 3}));

  Utils::Vector3i col = boost::qvm::col<0>(mat);
  BOOST_CHECK((col == Utils::Vector3i{1, 4, 7}));

  Utils::Matrix<int, 3, 3> diag_mat = boost::qvm::diag_mat(diagonal);
  BOOST_CHECK((diag_mat(0, 0) = 1));
  BOOST_CHECK((*(diag_mat.data()) == 1));
}

BOOST_AUTO_TEST_CASE(identity_matrix) {
  auto const id_mat = Utils::identity_mat<int, 2, 2>();
  BOOST_CHECK((id_mat(0, 0) == 1));
  BOOST_CHECK((id_mat(1, 1) == 1));
  BOOST_CHECK((id_mat(0, 1) == 0));
  BOOST_CHECK((id_mat(1, 0) == 0));
}

BOOST_AUTO_TEST_CASE(matrix_serialization) {
  Utils::Matrix<int, 2, 2> mat2{{8, 2}, {3, 4}};

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << mat2;

  Utils::Matrix<int, 2, 2> mat2_deserialized;
  boost::archive::text_iarchive in_ar(stream);
  in_ar >> mat2_deserialized;

  BOOST_CHECK((mat2_deserialized(0, 0) == 8));
  BOOST_CHECK((mat2_deserialized(1, 0) == 3));
  BOOST_CHECK((mat2_deserialized(0, 1) == 2));
  BOOST_CHECK((mat2_deserialized(1, 1) == 4));
}

BOOST_AUTO_TEST_CASE(type_deduction) {
  static_assert(
      std::is_same_v<boost::qvm::deduce_vec2<Utils::Matrix<int, 2, 2>,
                                             Utils::Vector<int, 2>, 2>::type,
                     Utils::Vector<int, 2>>);
  static_assert(
      std::is_same_v<boost::qvm::deduce_vec2<Utils::Matrix<int, 2, 2>,
                                             Utils::Vector<double, 2>, 2>::type,
                     Utils::Vector<double, 2>>);
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(matrix_matrix) {
  const Utils::Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Utils::Matrix<int, 2, 2> b{{5, 6}, {7, 8}};
  auto const res = a * b;
  BOOST_CHECK((res(0, 0) == 19));
  BOOST_CHECK((res(0, 1) == 22));
  BOOST_CHECK((res(1, 0) == 43));
  BOOST_CHECK((res(1, 1) == 50));
}

BOOST_AUTO_TEST_CASE(sub_add) {
  const Utils::Matrix<int, 2, 2> a{{1, 2}, {3, 4}};
  const Utils::Matrix<int, 2, 2> b{{5, 6}, {7, 8}};
  // unary -
  {
    auto const min = -a;
    BOOST_CHECK((min(0, 0) == -a(0, 0)));
    BOOST_CHECK((min(0, 1) == -a(0, 1)));
    BOOST_CHECK((min(1, 0) == -a(1, 0)));
    BOOST_CHECK((min(1, 1) == -a(1, 1)));
  }
  // binary -
  {
    auto const sub = b - a;
    BOOST_CHECK((sub(0, 0) == b(0, 0) - a(0, 0)));
    BOOST_CHECK((sub(0, 1) == b(0, 1) - a(0, 1)));
    BOOST_CHECK((sub(1, 0) == b(1, 0) - a(1, 0)));
    BOOST_CHECK((sub(1, 1) == b(1, 1) - a(1, 1)));
  }
  // binary +
  {
    auto const add = a + b;
    BOOST_CHECK((add(0, 0) == b(0, 0) + a(0, 0)));
    BOOST_CHECK((add(0, 1) == b(0, 1) + a(0, 1)));
    BOOST_CHECK((add(1, 0) == b(1, 0) + a(1, 0)));
    BOOST_CHECK((add(1, 1) == b(1, 1) + a(1, 1)));
  }
}

BOOST_AUTO_TEST_CASE(matrix_vector) {
  // trivial identity multiplication
  auto const id_mat = Utils::identity_mat<int, 2, 2>();
  auto const vec = Utils::Vector<int, 2>{2, 3};
  auto const res = id_mat * vec;
  BOOST_CHECK((res == vec));
  auto const res2 = id_mat.transposed() * vec;
  BOOST_CHECK((res2 == vec));

  // non-trivial multiplication
  auto const mat2 = Utils::Matrix<int, 2, 2>{{1, 2}, {3, 4}};
  auto const vec2 = Utils::Vector<int, 2>{7, 8};
  auto const res_vec = mat2 * vec2;
  BOOST_CHECK((res_vec[0] == 23));
  BOOST_CHECK((res_vec[1] == 53));
}

BOOST_AUTO_TEST_CASE(transposed) {
  auto const mat = Utils::Matrix<int, 2, 2>{{1, 2}, {3, 4}};
  auto const mat_T = mat.transposed();
  auto const expected = Utils::Matrix<int, 2, 2>{{1, 3}, {2, 4}};
  BOOST_CHECK((mat_T == expected));
}

BOOST_AUTO_TEST_CASE(inversed) {
  auto constexpr eps = std::numeric_limits<double>::epsilon();
  auto const mat = Utils::Matrix<double, 2, 2>{{1, 2}, {3, 4}};
  auto const mat_inv = mat.inversed();
  auto const expected = Utils::Matrix<double, 2, 2>{{-2, 1.0}, {1.5, -0.5}};
  BOOST_CHECK((mat_inv == expected));
  BOOST_CHECK_SMALL((mat_inv(0, 0) - (-2.0)), eps);
  BOOST_CHECK_SMALL((mat_inv(0, 1) - 1.0), eps);
  BOOST_CHECK_SMALL((mat_inv(1, 0) - 1.5), eps);
  BOOST_CHECK_SMALL((mat_inv(1, 1) - (-0.5)), eps);
}

BOOST_AUTO_TEST_CASE(shape) {
  auto const mat = Utils::Matrix<int, 2, 3>{{1, 2, 3}, {4, 5, 6}};
  BOOST_CHECK((mat.shape().first == 2));
  BOOST_CHECK((mat.shape().second == 3));
  BOOST_CHECK((mat.transposed().shape().first == 3));
  BOOST_CHECK((mat.transposed().shape().second == 2));
}

BOOST_AUTO_TEST_CASE(tensor_product) {
  // outer product of two vectors with same size
  auto const v1 = Utils::Vector<int, 2>{1, 2};
  auto const v2 = Utils::Vector<int, 2>{3, 4};
  Utils::Matrix<int, 2, 2> mat1 =
      boost::qvm::col_mat(v1) * boost::qvm::row_mat(v2);
  auto const expected = Utils::Matrix<int, 2, 2>{{3, 4}, {6, 8}};
  BOOST_CHECK((mat1 == expected));

  // outer product of two vectors with *different* size
  auto const v3 = Utils::Vector3i{5, 6, 7};
  Utils::Matrix<int, 2, 3> mat2 =
      boost::qvm::col_mat(v1) * boost::qvm::row_mat(v3);
  auto const expected2 = Utils::Matrix<int, 2, 3>{{5, 6, 7}, {10, 12, 14}};
  BOOST_CHECK((mat2 == expected2));
}

BOOST_AUTO_TEST_CASE(diagonal_matrix) {
  auto const v = Utils::Vector<int, 2>{1, 2};
  auto const m = Utils::diagonal_mat<int, 2, 2>(v);
  auto const expected = Utils::Matrix<int, 2, 2>{{1, 0}, {0, 2}};
  BOOST_CHECK((m == expected));
}

BOOST_AUTO_TEST_CASE(diagonal) {
  auto const v = Utils::Vector<int, 2>{1, 2};
  auto const m = Utils::diagonal_mat<int, 2, 2>(v);
  auto const diag = m.diagonal();
  BOOST_CHECK((v == diag));
}

BOOST_AUTO_TEST_CASE(trace) {
  auto const v = Utils::Vector<int, 2>{1, 2};
  auto const m = Utils::diagonal_mat<int, 2, 2>(v);
  auto const tr = m.trace();
  BOOST_CHECK((tr == 3));
}
