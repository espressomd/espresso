/*
 * Copyright (C) 2017-2026 The ESPResSo project
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
#define BOOST_TEST_MODULE VectorReference test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <array>
#include <cstddef>

// simulate a component-major (LayoutLeft) column with 4 rows:
// memory is [x0 x1 x2 x3 | y0 y1 y2 y3 | z0 z1 z2 z3], stride 4
struct ColumnFixture {
  std::array<double, 12> storage{};
  static constexpr std::size_t stride = 4;
  VectorReference row(std::size_t i) { return {storage.data() + i, stride}; }
};

BOOST_FIXTURE_TEST_CASE(assignment_writes_through, ColumnFixture) {
  row(1) = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(storage[1], 1.); // x1
  BOOST_CHECK_EQUAL(storage[5], 2.); // y1
  BOOST_CHECK_EQUAL(storage[9], 3.); // z1
}

BOOST_FIXTURE_TEST_CASE(conversion_reads_components, ColumnFixture) {
  storage[2] = 4.;  // x2
  storage[6] = 5.;  // y2
  storage[10] = 6.; // z2
  Utils::Vector3d const value = row(2);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 5.);
  BOOST_CHECK_EQUAL(value[2], 6.);
}

BOOST_FIXTURE_TEST_CASE(compound_operators, ColumnFixture) {
  row(0) = Utils::Vector3d{1., 1., 1.};
  row(0) += Utils::Vector3d{1., 2., 3.};
  row(0) -= Utils::Vector3d{0., 1., 0.};
  row(0) *= 2.;
  Utils::Vector3d const value = row(0);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 4.);
  BOOST_CHECK_EQUAL(value[2], 8.);
}

BOOST_FIXTURE_TEST_CASE(subscript_and_norms, ColumnFixture) {
  row(3)[0] = 3.;
  row(3)[1] = 0.;
  row(3)[2] = 4.;
  BOOST_CHECK_EQUAL(row(3)[2], 4.);
  BOOST_CHECK_EQUAL(row(3).norm2(), 25.);
  BOOST_CHECK_EQUAL(row(3).norm(), 5.);
}

BOOST_FIXTURE_TEST_CASE(proxy_copies_alias_the_same_row, ColumnFixture) {
  auto reference = row(1);
  auto alias = reference; // copying the proxy copies the reference, not data
  alias = Utils::Vector3d{7., 8., 9.};
  BOOST_CHECK_EQUAL(storage[1], 7.);
}

BOOST_FIXTURE_TEST_CASE(copy_assignment_copies_values, ColumnFixture) {
  row(0) = Utils::Vector3d{1., 2., 3.};
  row(1) = Utils::Vector3d{4., 5., 6.};
  // Copy-assignment must write VALUES through, not rebind: storage of row 0
  // changes to row 1's values, and row 1 itself is unchanged.
  row(0) = row(1);
  BOOST_CHECK_EQUAL(storage[0], 4.); // x0 <- x1
  BOOST_CHECK_EQUAL(storage[4], 5.); // y0 <- y1
  BOOST_CHECK_EQUAL(storage[8], 6.); // z0 <- z1
  // row 1 unchanged
  BOOST_CHECK_EQUAL(storage[1], 4.); // x1
  BOOST_CHECK_EQUAL(storage[5], 5.); // y1
  BOOST_CHECK_EQUAL(storage[9], 6.); // z1
  // writing through row(0) afterwards must not touch row 1 (proved distinct)
  row(0) = Utils::Vector3d{0., 0., 0.};
  BOOST_CHECK_EQUAL(storage[1], 4.); // x1 still row 1's value
}

BOOST_FIXTURE_TEST_CASE(chained_copy_assignment, ColumnFixture) {
  row(0) = Utils::Vector3d{9., 9., 9.};
  row(1) = Utils::Vector3d{9., 9., 9.};
  // Chained assignment must propagate the assigned values to both rows.
  row(0) = row(1) = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(storage[1], 1.); // x1
  BOOST_CHECK_EQUAL(storage[5], 2.); // y1
  BOOST_CHECK_EQUAL(storage[9], 3.); // z1
  BOOST_CHECK_EQUAL(storage[0], 1.); // x0
  BOOST_CHECK_EQUAL(storage[4], 2.); // y0
  BOOST_CHECK_EQUAL(storage[8], 3.); // z0
}

// ---------------------------------------------------------------------------
// IntegerVectorReference tests (BasicVectorReference<int>)
// ---------------------------------------------------------------------------

// simulate a component-major int column with 4 rows:
// memory is [x0 x1 x2 x3 | y0 y1 y2 y3 | z0 z1 z2 z3], stride 4
struct IntColumnFixture {
  std::array<int, 12> storage{};
  static constexpr std::size_t stride = 4;
  IntegerVectorReference row(std::size_t i) {
    return {storage.data() + i, stride};
  }
};

BOOST_FIXTURE_TEST_CASE(int_assignment_writes_through, IntColumnFixture) {
  row(1) = Utils::Vector<int, 3>{10, 20, 30};
  BOOST_CHECK_EQUAL(storage[1], 10); // x1
  BOOST_CHECK_EQUAL(storage[5], 20); // y1
  BOOST_CHECK_EQUAL(storage[9], 30); // z1
}

BOOST_FIXTURE_TEST_CASE(int_conversion_reads_components, IntColumnFixture) {
  storage[2] = 4;  // x2
  storage[6] = 5;  // y2
  storage[10] = 6; // z2
  Utils::Vector<int, 3> const value = row(2);
  BOOST_CHECK_EQUAL(value[0], 4);
  BOOST_CHECK_EQUAL(value[1], 5);
  BOOST_CHECK_EQUAL(value[2], 6);
}

BOOST_FIXTURE_TEST_CASE(int_subscript_read_write, IntColumnFixture) {
  row(3)[0] = 7;
  row(3)[1] = 8;
  row(3)[2] = 9;
  BOOST_CHECK_EQUAL(row(3)[0], 7);
  BOOST_CHECK_EQUAL(row(3)[1], 8);
  BOOST_CHECK_EQUAL(row(3)[2], 9);
  // verify storage layout
  BOOST_CHECK_EQUAL(storage[3], 7);  // x3
  BOOST_CHECK_EQUAL(storage[7], 8);  // y3
  BOOST_CHECK_EQUAL(storage[11], 9); // z3
}

// ---------------------------------------------------------------------------
// QuaternionReference tests
// ---------------------------------------------------------------------------

// simulate a component-major double column for quaternions with 4 rows:
// memory is [w0 w1 w2 w3 | x0 x1 x2 x3 | y0 y1 y2 y3 | z0 z1 z2 z3],
// stride 4
struct QuatColumnFixture {
  std::array<double, 16> storage{};
  static constexpr std::size_t stride = 4;
  QuaternionReference row(std::size_t i) {
    return {storage.data() + i, stride};
  }
};

BOOST_FIXTURE_TEST_CASE(quat_assignment_writes_through, QuatColumnFixture) {
  Utils::Quaternion<double> q;
  q[0] = 1.0;
  q[1] = 2.0;
  q[2] = 3.0;
  q[3] = 4.0;
  row(1) = q;
  BOOST_CHECK_EQUAL(storage[1], 1.0);  // q[0] at row 1
  BOOST_CHECK_EQUAL(storage[5], 2.0);  // q[1] at row 1
  BOOST_CHECK_EQUAL(storage[9], 3.0);  // q[2] at row 1
  BOOST_CHECK_EQUAL(storage[13], 4.0); // q[3] at row 1
}

BOOST_FIXTURE_TEST_CASE(quat_subscript_read_write, QuatColumnFixture) {
  row(2)[0] = 0.5;
  row(2)[1] = 0.5;
  row(2)[2] = 0.5;
  row(2)[3] = 0.5;
  BOOST_CHECK_EQUAL(row(2)[0], 0.5);
  BOOST_CHECK_EQUAL(row(2)[1], 0.5);
  BOOST_CHECK_EQUAL(row(2)[2], 0.5);
  BOOST_CHECK_EQUAL(row(2)[3], 0.5);
  // verify storage layout
  BOOST_CHECK_EQUAL(storage[2], 0.5);  // component 0 of row 2
  BOOST_CHECK_EQUAL(storage[6], 0.5);  // component 1 of row 2
  BOOST_CHECK_EQUAL(storage[10], 0.5); // component 2 of row 2
  BOOST_CHECK_EQUAL(storage[14], 0.5); // component 3 of row 2
}

BOOST_FIXTURE_TEST_CASE(quat_plus_equals, QuatColumnFixture) {
  row(0)[0] = 1.0;
  row(0)[1] = 2.0;
  row(0)[2] = 3.0;
  row(0)[3] = 4.0;
  Utils::Quaternion<double> delta;
  delta[0] = 0.1;
  delta[1] = 0.2;
  delta[2] = 0.3;
  delta[3] = 0.4;
  row(0) += delta;
  BOOST_CHECK_CLOSE(row(0)[0], 1.1, 1e-9);
  BOOST_CHECK_CLOSE(row(0)[1], 2.2, 1e-9);
  BOOST_CHECK_CLOSE(row(0)[2], 3.3, 1e-9);
  BOOST_CHECK_CLOSE(row(0)[3], 4.4, 1e-9);
}

BOOST_FIXTURE_TEST_CASE(quat_divide_equals, QuatColumnFixture) {
  row(0)[0] = 2.0;
  row(0)[1] = 4.0;
  row(0)[2] = 6.0;
  row(0)[3] = 8.0;
  row(0) /= 2.0;
  BOOST_CHECK_EQUAL(row(0)[0], 1.0);
  BOOST_CHECK_EQUAL(row(0)[1], 2.0);
  BOOST_CHECK_EQUAL(row(0)[2], 3.0);
  BOOST_CHECK_EQUAL(row(0)[3], 4.0);
}

BOOST_FIXTURE_TEST_CASE(quat_conversion_round_trip, QuatColumnFixture) {
  Utils::Quaternion<double> q;
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
  row(0) = q;
  Utils::Quaternion<double> const result = row(0);
  BOOST_CHECK_EQUAL(result[0], 1.0);
  BOOST_CHECK_EQUAL(result[1], 0.0);
  BOOST_CHECK_EQUAL(result[2], 0.0);
  BOOST_CHECK_EQUAL(result[3], 0.0);
}

BOOST_FIXTURE_TEST_CASE(quat_copy_assignment_writes_through,
                        QuatColumnFixture) {
  Utils::Quaternion<double> q0;
  q0[0] = 1.0;
  q0[1] = 2.0;
  q0[2] = 3.0;
  q0[3] = 4.0;
  Utils::Quaternion<double> q1;
  q1[0] = 5.0;
  q1[1] = 6.0;
  q1[2] = 7.0;
  q1[3] = 8.0;
  row(0) = q0;
  row(1) = q1;
  // copy-assignment must write VALUES through, not rebind
  row(0) = row(1);
  // row 0 storage should now hold row 1's values
  BOOST_CHECK_EQUAL(storage[0], 5.0);  // component 0 of row 0
  BOOST_CHECK_EQUAL(storage[4], 6.0);  // component 1 of row 0
  BOOST_CHECK_EQUAL(storage[8], 7.0);  // component 2 of row 0
  BOOST_CHECK_EQUAL(storage[12], 8.0); // component 3 of row 0
  // row 1 unchanged
  BOOST_CHECK_EQUAL(storage[1], 5.0);  // component 0 of row 1
  BOOST_CHECK_EQUAL(storage[5], 6.0);  // component 1 of row 1
  BOOST_CHECK_EQUAL(storage[9], 7.0);  // component 2 of row 1
  BOOST_CHECK_EQUAL(storage[13], 8.0); // component 3 of row 1
  // writing through row(0) afterwards must not touch row 1
  row(0) = Utils::Quaternion<double>{};
  BOOST_CHECK_EQUAL(storage[1], 5.0); // component 0 of row 1 still unchanged
}

BOOST_FIXTURE_TEST_CASE(quat_norm_and_norm2, QuatColumnFixture) {
  // q = (1, 0, 0, 0) has norm 1 and norm2 1
  row(0)[0] = 1.0;
  row(0)[1] = 0.0;
  row(0)[2] = 0.0;
  row(0)[3] = 0.0;
  BOOST_CHECK_EQUAL(row(0).norm2(), 1.0);
  BOOST_CHECK_EQUAL(row(0).norm(), 1.0);
  // q = (0, 2, 0, 0) has norm 2 and norm2 4
  row(1)[0] = 0.0;
  row(1)[1] = 2.0;
  row(1)[2] = 0.0;
  row(1)[3] = 0.0;
  BOOST_CHECK_EQUAL(row(1).norm2(), 4.0);
  BOOST_CHECK_EQUAL(row(1).norm(), 2.0);
}
