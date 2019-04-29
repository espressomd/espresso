/*
  Copyright (C) 2017-2018 The ESPResSo project

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
#define BOOST_TEST_MODULE all_compare test
#define BOOST_TEST_DYN_LINK
#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/mpi/gatherv.hpp"
using Utils::Mpi::gatherv;

#include <string>

namespace mpi = boost::mpi;

/*
 * Check that implementation behaves
 * like MPI_Gatherv with an mpi datatype.
 * To test this we gather the rank from
 * every rank to one rank and then check
 * that the value was written to the
 * correct position in the output array.
 */
BOOST_AUTO_TEST_CASE(mpi_type) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const root = world.size() - 1;

  /* out-of-place */
  {
    if (rank == root) {
      std::vector<int> out(size, -1);
      std::vector<int> sizes(size, 1);

      gatherv(world, &rank, 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(i, out.at(i));
      }
    } else {
      Utils::Mpi::gatherv(world, &rank, 1, root);
    }
  }

  /* in-place */
  {
    if (rank == root) {
      std::vector<int> out(size, -1);
      out[rank] = rank;
      std::vector<int> sizes(size, 1);

      gatherv(world, out.data(), 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(i, out.at(i));
      }
    } else {
      Utils::Mpi::gatherv(world, &rank, 1, root);
    }
  }
}

/*
 * Check that implementation behaves
 * like MPI_Gatherv with an non-mpi datatype.
 * To test this we gather a string containing the rank from
 * every rank to one rank and then check
 * that the value was written to the
 * correct position in the output array.
 */
BOOST_AUTO_TEST_CASE(non_mpi_type) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const root = world.size() - 1;
  auto const in = std::to_string(rank);

  /* out-of-place */
  {
    if (rank == root) {
      std::vector<std::string> out(size);
      std::vector<int> sizes(size, 1);

      gatherv(world, &in, 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(std::to_string(i), out.at(i));
      }
    } else {
      Utils::Mpi::gatherv(world, &in, 1, root);
    }
  }

  /* in-place */
  {
    if (rank == root) {
      std::vector<std::string> out(size);
      out[rank] = in;
      std::vector<int> sizes(size, 1);

      gatherv(world, out.data(), 1, out.data(), sizes.data(), root);

      for (int i = 0; i < size; i++) {
        BOOST_CHECK_EQUAL(std::to_string(i), out.at(i));
      }
    } else {
      Utils::Mpi::gatherv(world, &in, 1, root);
    }
  }
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
