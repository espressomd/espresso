/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

/* Unit tests for the ErrorHandling::RuntimeErrorCollector class. */

#include <algorithm>
#include <memory>

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE RuntimeError test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "RuntimeErrorCollector.hpp"

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

namespace Testing {

void reduce_and_check(const boost::mpi::communicator &comm, bool local_value) {
  if (comm.rank() == 0) {
    bool total = true;
    boost::mpi::reduce(comm, local_value, total, std::logical_and<bool>(), 0);
    BOOST_CHECK(total);
  } else {
    boost::mpi::reduce(comm, local_value, std::logical_and<bool>(), 0);
  }
}
} // namespace Testing

using ErrorHandling::RuntimeError;
using ErrorHandling::RuntimeErrorCollector;

BOOST_AUTO_TEST_CASE(count) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  Testing::reduce_and_check(world, rec.count() == 0);

  /* MPI guarantees that size >= 1 and rank 0 exists. */
  if (world.rank() == (world.size() - 1)) {
    rec.error("Test_error", "Test_functions", "Test_file", 42);
  }

  Testing::reduce_and_check(world, rec.count() == 1);

  rec.warning("Test_error", "Test_functions", "Test_file", 42);

  Testing::reduce_and_check(world, rec.count() == world.size() + 1);

  /* There should now be one error and world.size() warnings */
  Testing::reduce_and_check(world,
                            rec.count(RuntimeError::ErrorLevel::ERROR) == 1);
  /* All messages are at least WARNING or higher. */
  {
    /* Beware of the execution order */
    int total = rec.count();
    Testing::reduce_and_check(
        world, rec.count(RuntimeError::ErrorLevel::WARNING) == total);
  }
}

/*
 * Check the message gathering. Every node generates a runtime error
 * and a warning. Then we gather the messages
 * on the master and check if we got the correct messages. Then we
 * check the post-condition count() == 0.
 */
BOOST_AUTO_TEST_CASE(gather) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  rec.error("Test_error", "Test_functions", "Test_file", world.rank());
  rec.warning("Test_error", "Test_functions", "Test_file", world.rank());

  if (world.rank() == 0) {
    /* Gathered error messages */
    auto results = rec.gather();
    /* Track how many messages we have seen from which node. */
    std::vector<int> present(world.size());

    for (const auto &err : results) {
      present[err.who()]++;
    }

    /* Check if we got 2 messages from every node. */
    BOOST_CHECK(std::all_of(present.begin(), present.end(),
                            [](int i) { return i == 2; }));
    /* Count warnings, should be world.rank() many */
    BOOST_CHECK(std::count_if(
                    results.begin(), results.end(), [](const RuntimeError &e) {
                      return e.level() == RuntimeError::ErrorLevel::WARNING;
                    }) == world.size());
    /* Count errors, should be world.rank() many */
    BOOST_CHECK(std::count_if(
                    results.begin(), results.end(), [](const RuntimeError &e) {
                      return e.level() == RuntimeError::ErrorLevel::ERROR;
                    }) == world.size());
  } else {
    rec.gatherSlave();
  }

  Testing::reduce_and_check(world, rec.count() == 0);
}
