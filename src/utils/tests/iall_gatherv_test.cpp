/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Utils::Mpi::iall_gatherv test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/mpi/iall_gatherv.hpp>

#include <boost/mpi.hpp>

#include <algorithm>
#include <string>
#include <vector>

/*
 * Check that implementation behaves like @c MPI_Iallgatherv.
 * To test this, we gather one value from
 * every rank to all ranks and then check
 * that the value was written to the
 * correct position in the output array.
 */
BOOST_AUTO_TEST_CASE(one_element) {
  boost::mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();

  /* out-of-place */
  {
    std::vector<int> out(size, -1);
    std::vector<int> sizes(size, 1);
    auto reqs =
        Utils::Mpi::iall_gatherv(world, &rank, 1, out.data(), sizes.data());
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(out.at(i), i);
    }
  }

  /* in-place */
  {
    std::vector<int> out(size, -1);
    std::vector<int> sizes(size, 1);
    out[rank] = rank;
    auto reqs = Utils::Mpi::iall_gatherv(world, out.data(), 1, out.data(),
                                         sizes.data());
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(out.at(i), i);
    }
  }
}

/*
 * Check that implementation behaves like @c MPI_Iallgatherv.
 * To test this, we gather values from
 * every rank to all ranks and then check
 * that the values were written to the
 * correct positions in the output array.
 */
BOOST_AUTO_TEST_CASE(multiple_elements) {
  boost::mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();

  std::vector<int> expected_values;
  for (int i = 0; i < size; ++i) {
    for (auto j = 0; j < i + 1; ++j) {
      expected_values.push_back(i);
    }
  }

  /* out-of-place */
  {
    std::vector<int> in(rank + 1, rank);
    std::vector<int> out(size * (size + 1) / 2, -1);
    std::vector<int> sizes(size, -1);
    std::iota(sizes.begin(), sizes.end(), 1);
    auto reqs = Utils::Mpi::iall_gatherv(world, in.data(), rank + 1, out.data(),
                                         sizes.data());
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    for (int i = 0; i < expected_values.size(); i++) {
      BOOST_CHECK_EQUAL(out.at(i), expected_values.at(i));
    }
  }

  /* in-place */
  {
    std::vector<int> out(size * (size + 1) / 2, -1);
    std::vector<int> sizes(size, -1);
    std::iota(sizes.begin(), sizes.end(), 1);
    for (int i = 0; i < rank + 1; ++i) {
      out[rank * (rank + 1) / 2 + i] = rank;
    }
    auto reqs = Utils::Mpi::iall_gatherv(world, out.data(), rank + 1,
                                         out.data(), sizes.data());
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    for (int i = 0; i < expected_values.size(); i++) {
      BOOST_CHECK_EQUAL(out.at(i), expected_values.at(i));
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
