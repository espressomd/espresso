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
#define BOOST_TEST_MODULE ScriptInterface::ParallelExceptionHandler test
#define BOOST_TEST_DYN_LINK

/* Guard against a GCC 12 diagnostic for Boost versions <= 1.79.
 * More details in https://github.com/boostorg/function/issues/42
 */
#include <boost/serialization/version.hpp>
#if BOOST_VERSION <= 107900 and defined(BOOST_GCC) and (BOOST_GCC >= 120000)
#define BOOST_HAS_GCC_12_UNINITIALIZED_DIAGNOSTIC
#endif
#if defined(BOOST_HAS_GCC_12_UNINITIALIZED_DIAGNOSTIC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <boost/test/unit_test.hpp>
#if defined(BOOST_HAS_GCC_12_UNINITIALIZED_DIAGNOSTIC)
#pragma GCC diagnostic pop
#endif

#include "script_interface/Exception.hpp"
#include "script_interface/ParallelExceptionHandler.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/errorhandling.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace utf = boost::unit_test;

static std::weak_ptr<boost::mpi::environment> mpi_env;

namespace Testing {
struct Error : public std::exception {};
struct Warning : public std::exception {};
} // namespace Testing

/** Decorator to skip tests running on only 1 MPI rank. */
struct if_parallel_test {
  boost::test_tools::assertion_result operator()(utf::test_unit_id) {
    boost::mpi::communicator world;
    return world.size() >= 2;
  }
};

BOOST_TEST_DECORATOR(*utf::precondition(if_parallel_test()))
BOOST_AUTO_TEST_CASE(parallel_exceptions) {
  boost::mpi::communicator world;
  auto callbacks =
      std::make_shared<Communication::MpiCallbacks>(world, ::mpi_env.lock());
  ErrorHandling::init_error_handling(callbacks);
  auto handler = ScriptInterface::ParallelExceptionHandler{world};

  {
    // exception on main rank -> re-throw on main rank
    bool rethrown = false;
    bool converted = false;
    try {
      handler.parallel_try_catch<Testing::Error>(
          []() { throw Testing::Error(); });
    } catch (Testing::Error const &err) {
      rethrown = true;
    } catch (ScriptInterface::Exception const &err) {
      converted = true;
    }
    if (world.rank() == 0) {
      BOOST_CHECK(rethrown);
    } else {
      BOOST_CHECK(converted);
    }
  }
  {
    // exception of an unknown type: not caught
    bool unhandled = false;
    try {
      handler.parallel_try_catch<Testing::Error>(
          []() { throw Testing::Warning(); });
    } catch (Testing::Warning const &err) {
      unhandled = true;
    }
    BOOST_CHECK(unhandled);
  }
  {
    // exception on worker rank -> communicate to main rank
    bool communicated = false;
    bool converted = false;
    try {
      handler.parallel_try_catch<Testing::Error>([&world]() {
        runtimeErrorMsg() << "harmless message";
        if (world.rank() != 0) {
          throw Testing::Error();
        }
      });
    } catch (std::runtime_error const &err) {
      communicated = true;
      if (world.rank() == 0) {
        BOOST_CHECK_EQUAL(err.what(),
                          "an error occurred on one or more MPI ranks:\n  rank "
                          "0: \n  rank 1: std::exception");
      }
    } catch (ScriptInterface::Exception const &err) {
      converted = true;
    }
    if (world.rank() == 0) {
      BOOST_CHECK(communicated);
    } else {
      BOOST_CHECK(converted);
    }
    // runtime error messages are printed to stderr and cleared
    BOOST_CHECK_EQUAL(check_runtime_errors_local(), 0);
  }
  if (world.rank() != 0) {
    callbacks->loop();
  }
}

int main(int argc, char **argv) {
  auto const mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  ::mpi_env = mpi_env;

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
