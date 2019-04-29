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
#define BOOST_TEST_MODULE all_gather test
#define BOOST_TEST_DYN_LINK
#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/mpi/all_gatherv.hpp"
using Utils::Mpi::all_gatherv;

#include <string>

namespace mpi = boost::mpi;

BOOST_AUTO_TEST_CASE(mpi_type) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();

  /* out-of-place */
  {
    std::vector<int> out(size, -1);
    std::vector<int> sizes(size, 1);

    all_gatherv(world, &rank, 1, out.data(), sizes.data());

    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(i, out.at(i));
    }
  }

  /* in-place */
  {
    std::vector<int> out(size, -1);
    out[rank] = rank;
    std::vector<int> sizes(size, 1);

    all_gatherv(world, out.data(), 1, out.data(), sizes.data());

    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(i, out.at(i));
    }
  }
}

BOOST_AUTO_TEST_CASE(non_mpi_type) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const in = std::to_string(rank);

  /* out-of-place */
  {
    std::vector<std::string> out(size);
    std::vector<int> sizes(size, 1);

    all_gatherv(world, &in, 1, out.data(), sizes.data());

    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(std::to_string(i), out.at(i));
    }
  }

  /* in-place */
  {
    std::vector<std::string> out(size);
    out[rank] = in;
    std::vector<int> sizes(size, 1);

    all_gatherv(world, out.data(), 1, out.data(), sizes.data());

    for (int i = 0; i < size; i++) {
      BOOST_CHECK_EQUAL(std::to_string(i), out.at(i));
    }
  }
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
