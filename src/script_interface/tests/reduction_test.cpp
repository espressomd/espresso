/*
 * Copyright (C) 2022 The ESPResSo project
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
#define BOOST_TEST_MODULE MPI reduction algorithms
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/communication.hpp"

#include <boost/mpi.hpp>
#include <boost/optional.hpp>

BOOST_AUTO_TEST_CASE(reduce_sum) {
  boost::mpi::communicator comm;
  auto const sum = ScriptInterface::mpi_reduce_sum(comm, comm.rank() + 1);

  if (comm.rank() == 0) {
    auto const n = comm.size();
    BOOST_CHECK_EQUAL(sum, (n * (n + 1)) / 2);
  } else {
    BOOST_CHECK_EQUAL(sum, 0);
  }
}

BOOST_AUTO_TEST_CASE(reduce_optional) {
  boost::mpi::communicator comm;

  for (int rank = 0; rank < comm.size(); ++rank) {
    boost::optional<int> maybe;
    if (comm.rank() == rank) {
      maybe = 42;
    }
    auto const sum = ScriptInterface::mpi_reduce_optional(comm, maybe);

    if (comm.rank() == 0) {
      BOOST_CHECK_EQUAL(sum, 42);
    } else {
      BOOST_CHECK_EQUAL(sum, 0);
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
