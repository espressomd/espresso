/*
 * Copyright (C) 2026 The ESPResSo project
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

/* Unit tests for the ErrorHandling::RuntimeErrorStream class. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE RuntimeErrorStream test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "error_handling/RuntimeError.hpp"
#include "error_handling/RuntimeErrorCollector.hpp"
#include "error_handling/RuntimeErrorStream.hpp"

#include <boost/mpi.hpp>

#include <cstddef>
#include <cstring>
#include <new>
#include <string>

using ErrorHandling::RuntimeError;
using ErrorHandling::RuntimeErrorCollector;
using ErrorHandling::RuntimeErrorStream;
using ErrorLevel = ErrorHandling::RuntimeError::ErrorLevel;

/*
 * Regression test for the RuntimeErrorStream copy constructor.
 *
 * The user-defined copy constructor must copy every data member.
 *
 * Test mechanics:
 *  - The copy is placement-new constructed into a raw stack buffer. Before
 *    construction we deterministically write a *valid but non-ERROR* byte
 *    pattern (WARNING == 0) over the whole buffer.
 *  - The copy constructor copies all data members, and thus the error
 *    collector must contain two exact copies of the error message
 */
BOOST_AUTO_TEST_CASE(copy_ctor_preserves_level) {
  boost::mpi::communicator world;

  RuntimeErrorCollector rec(world);

  std::string const unique_file = "unique_file_marker";
  std::string const source_text = "source_text";

  if (world.rank() == world.size() - 1) {
    // Raw aligned storage, deterministically pre-filled with bytes that
    // encode the valid ErrorLevel WARNING (== 0). On a broken copy
    // constructor m_level is never written, so it reads back as WARNING.
    alignas(
        RuntimeErrorStream) unsigned char storage[sizeof(RuntimeErrorStream)];
    static_assert(static_cast<int>(ErrorLevel::WARNING) == 0,
                  "test assumes WARNING maps to all-zero bytes");
    std::memset(storage, 0, sizeof(storage));

    // Source stream: ERROR level, unique file, non-empty marker text.
    RuntimeErrorStream source(rec, ErrorLevel::ERROR, unique_file, 42,
                              "Test_function");
    source << source_text;

    // Copy-construct into the pre-filled storage. This invokes the copy
    // constructor on a genuine lvalue (never elided).
    auto *copy = new (storage) RuntimeErrorStream(source);

    // Manually invoke the destructor: it reads m_level and records a message
    copy->~RuntimeErrorStream();
    // The `source` destructor is invoked at end of scope and records a message
  }

  world.barrier();

  if (world.rank() == 0) {
    auto const results = rec.gather();
    BOOST_REQUIRE_EQUAL(results.size(), 2ul);
    for (auto const &err : results) {
      BOOST_CHECK_EQUAL(static_cast<int>(err.level()),
                        static_cast<int>(ErrorLevel::ERROR));
      BOOST_CHECK_EQUAL(err.what(), source_text);
      BOOST_CHECK_EQUAL(err.who(), world.size() - 1);
      BOOST_CHECK_EQUAL(err.function(), "Test_function");
      BOOST_CHECK_EQUAL(err.file(), unique_file);
      BOOST_CHECK_EQUAL(err.line(), 42);
    }
  } else {
    rec.gather_local();
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
