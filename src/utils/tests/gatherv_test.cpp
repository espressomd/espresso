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
#define BOOST_TEST_MODULE Utils::Mpi::gatherv test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/mpi/gatherv.hpp>

#include <boost/mpi.hpp>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

struct identity {
  template <class T> constexpr T &&operator()(T &&t) const noexcept {
    return std::forward<T>(t);
  }
};

template <typename T, class F> void check(T default_value, F conversion) {
  boost::mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const root = world.size() - 1;
  auto const in = conversion(rank);

  /* out-of-place */
  {
    if (rank == root) {
      std::vector<T> out(size, default_value);
      std::vector<int> sizes(size, 1);

      Utils::Mpi::gatherv(world, &in, 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(out.at(i), conversion(i));
      }
    } else if (rank % 2 == 0) {
      Utils::Mpi::gatherv(world, &in, 1, static_cast<T *>(nullptr), nullptr,
                          root);
    } else {
      Utils::Mpi::gatherv(world, &in, 1, root);
    }
  }

  /* in-place */
  {
    if (rank == root) {
      std::vector<T> out(size, default_value);
      out[rank] = in;
      std::vector<int> sizes(size, 1);

      Utils::Mpi::gatherv(world, out.data(), 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(out.at(i), conversion(i));
      }
    } else if (rank % 2 == 0) {
      Utils::Mpi::gatherv(world, &in, 1, static_cast<T *>(nullptr), nullptr,
                          root);
    } else {
      Utils::Mpi::gatherv(world, &in, 1, root);
    }
  }
}

/*
 * Check that implementation behaves
 * like @c MPI_Gatherv with an mpi datatype.
 * To test this, we gather the rank from
 * every rank to one rank and then check
 * that the value was written to the
 * correct position in the output array.
 */
BOOST_AUTO_TEST_CASE(mpi_type) { check(-1, identity{}); }

/*
 * Check that implementation behaves
 * like @c MPI_Gatherv with a non-mpi datatype.
 * To test this, we gather a string containing the rank from
 * every rank to one rank and then check
 * that the value was written to the
 * correct position in the output array.
 */
BOOST_AUTO_TEST_CASE(non_mpi_type) {
  check(std::string{""}, static_cast<std::string (&)(int)>(std::to_string));
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
