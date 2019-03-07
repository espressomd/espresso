/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Utils::Array test
#define BOOST_TEST_DYN_LINK

#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/Array.hpp"
using Utils::Array;

#include <array>
#include <numeric>
#include <vector>

BOOST_AUTO_TEST_CASE(const_expr_ctor) {
  static_assert(4 == Array<int, 4>().size(), "");
  static_assert(4 == Array<int, 4>().max_size(), "");
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
  auto a = Array<int, 4>{{1, 2, 3, 4}};

  BOOST_CHECK(*(a.begin()) == 1);
  BOOST_CHECK(*(a.cbegin()) == 1);
  BOOST_CHECK(*(a.end() - 1) == 4);
  BOOST_CHECK(*(a.cend() - 1) == 4);
}

BOOST_AUTO_TEST_CASE(element_access) {
  auto a = Array<int, 5>{{5, 6, 7, 8, 9}};

  int c = 5;
  for (Array<int, 5>::size_type i = 0; i < a.size(); ++i) {
    BOOST_CHECK(a.at(i) == c);
    BOOST_CHECK(a[i] == c);
    ++c;
  }

  BOOST_CHECK_THROW(a.at(a.size()), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(fill) {
  Array<int, 10> a{};
  a.fill(10);
  for (Array<int, 10>::size_type i = 0; i < a.size(); ++i) {
    BOOST_CHECK(a[i] == 10);
  }
}

BOOST_AUTO_TEST_CASE(broadcast) {
  constexpr auto a = Array<int, 3>::broadcast(5);
  static_assert(a[0] == 5, "");
  static_assert(a[1] == 5, "");
  static_assert(a[2] == 5, "");
}

BOOST_AUTO_TEST_CASE(serialization) {
  boost::mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  Array<double, 3> a{1.0, 2.0, 3.0};
  if (size > 1) {
    if (rank == 0) {
      world.send(1, 42, a);
    } else if (rank == 1) {
      Array<double, 3> b{};
      world.recv(0, 42, b);
      for (Array<double, 3>::size_type i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(a[i], b[i]);
      }
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
