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
#define BOOST_TEST_MODULE sendrecv test
#define BOOST_TEST_DYN_LINK

#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include <string>

#include "utils/mpi/sendrecv.hpp"

using Utils::Mpi::isendrecv;
using Utils::Mpi::sendrecv;

namespace mpi = boost::mpi;

BOOST_AUTO_TEST_CASE(blocking) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const left = (rank - 1 + size) % size;
  auto const right = (rank + 1) % size;

  /* MPI datatype */
  {
    const int send = rank;
    int recv;

    sendrecv(world, left, 42, send, right, 42, recv);

    BOOST_CHECK_EQUAL(recv, right);
  }

  /* Non-MPI datatype */
  {
    const std::string send = std::to_string(rank);
    std::string recv;

    sendrecv(world, left, 42, send, right, 42, recv);

    BOOST_CHECK_EQUAL(recv, std::to_string(right));
  }
}

BOOST_AUTO_TEST_CASE(non_blocking) {
  mpi::communicator world;
  auto const rank = world.rank();
  auto const size = world.size();
  auto const left = (rank - 1 + size) % size;
  auto const right = (rank + 1) % size;

  const int send = rank;
  int recv;

  auto req = isendrecv(world, left, 42, send, right, 42, recv);
  mpi::wait_all(req.begin(), req.end());

  BOOST_CHECK_EQUAL(recv, right);
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
