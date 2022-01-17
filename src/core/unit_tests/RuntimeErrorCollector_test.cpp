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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE RuntimeError test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "RuntimeErrorCollector.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <vector>

using ErrorHandling::RuntimeError;
using ErrorHandling::RuntimeErrorCollector;
using ErrorLevel = ErrorHandling::RuntimeError::ErrorLevel;

BOOST_AUTO_TEST_CASE(count) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  BOOST_REQUIRE_EQUAL(rec.count(), 0);
  world.barrier();

  /* MPI guarantees that size >= 1 and rank 0 exists. */
  auto const rank_of_error = world.size() - 1;
  if (world.rank() == rank_of_error) {
    rec.error("Test_error", "Test_functions", "Test_file", 42);
  }
  world.barrier();
  BOOST_REQUIRE_EQUAL(rec.count(), 1);
  world.barrier();
  rec.warning("Test_error", "Test_functions", "Test_file", 42);
  world.barrier();
  BOOST_REQUIRE_EQUAL(rec.count(), world.size() + 1);

  /* There should now be one error and world.size() warnings */
  auto const n_local_errors = static_cast<int>(world.rank() == rank_of_error);
  auto const n_local_warnings = n_local_errors + 1;
  BOOST_REQUIRE_EQUAL(rec.count(ErrorLevel::ERROR), n_local_errors);
  BOOST_REQUIRE_EQUAL(rec.count(ErrorLevel::WARNING), n_local_warnings);

  /* Clear list of errors */
  world.barrier();
  rec.clear();
  world.barrier();
  BOOST_REQUIRE_EQUAL(rec.count(), 0);
}

/*
 * Check the message gathering. Every node generates a runtime error
 * and a warning. Then we gather the messages
 * on the head node and check if we got the correct messages. Then we
 * check the post-condition count() == 0.
 */
BOOST_AUTO_TEST_CASE(gather) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  rec.error("Test_error", "Test_functions", "Test_file", world.rank());
  rec.warning("Test_error", "Test_functions", "Test_file", world.rank());
  world.barrier();

  if (world.rank() == 0) {
    /* Gathered error messages */
    auto const results = rec.gather();
    /* Track how many messages we have seen from which node. */
    std::vector<int> present(world.size());

    for (const auto &err : results) {
      present[err.who()]++;
    }

    /* Check if we got 2 messages from every node. */
    BOOST_CHECK(std::all_of(present.begin(), present.end(),
                            [](int i) { return i == 2; }));
    /* Count warnings, should be world.rank() many */
    BOOST_CHECK(std::count_if(results.begin(), results.end(),
                              [](const RuntimeError &e) {
                                return e.level() == ErrorLevel::WARNING;
                              }) == world.size());
    /* Count errors, should be world.rank() many */
    BOOST_CHECK(std::count_if(results.begin(), results.end(),
                              [](const RuntimeError &e) {
                                return e.level() == ErrorLevel::ERROR;
                              }) == world.size());
  } else {
    rec.gather_local();
  }

  BOOST_REQUIRE_EQUAL(rec.count(), 0);
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
