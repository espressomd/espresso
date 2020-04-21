/*
 * Copyright (C) 2017-2019 The ESPResSo project
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
#define BOOST_TEST_MODULE all_compare test
#define BOOST_TEST_DYN_LINK
#include <boost/mpi.hpp>
#include <boost/test/unit_test.hpp>

#include "utils/mpi/all_compare.hpp"
using Utils::Mpi::all_compare;

namespace mpi = boost::mpi;

BOOST_AUTO_TEST_CASE(true_) {
  mpi::communicator world;

  BOOST_CHECK(all_compare(world, 42));
}

BOOST_AUTO_TEST_CASE(false_) {
  mpi::communicator world;

  BOOST_CHECK(all_compare(world, (world.rank() > 0) ? 42 : 41) ==
              (world.size() <= 1));
}

int main(int argc, char **argv) {
  mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
