/*
  Copyright (C) 2016 The ESPResSo project
    Max-Planck-Institute for Polymer Research, Theory Group

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

/** \file MpiCallbacks_test.cpp Unit tests for the MpiCallbacks class.
 *
*/

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE MpiCallbacks test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../MpiCallbacks.hpp"

namespace Testing {
void reduce_and_check(const boost::mpi::communicator &comm, bool local_value) {
  if (comm.rank() == 0) {
    bool total;
    boost::mpi::reduce(comm, local_value, total, std::logical_and<bool>(), 0);
    BOOST_CHECK(total); // NOLINT
  } else {
    boost::mpi::reduce(comm, local_value, std::logical_and<bool>(), 0);
  }
}
}

bool static_callback_called = false;
void static_callback(int a, int b) {
  static_callback_called = (a == 5) && (b == 13);
}

using boost::mpi::communicator;
using Communication::MpiCallbacks;

/**
 * Check is the mpi loop can be aborted
 * by abort_loop()
 */
BOOST_AUTO_TEST_CASE(loop_exit) {
  communicator world;
  MpiCallbacks callbacks(world, /* abort_on_exit */ false);

  if (world.rank() == 0) {
    callbacks.abort_loop();
  } else {
    callbacks.loop();
  }

  Testing::reduce_and_check(world, true);
}

/**
 * Check if adding and calling a static callback
 * works as expected.
 */
BOOST_AUTO_TEST_CASE(add_static_callback) {
  communicator world;
  MpiCallbacks callbacks(world, /* abort_on_exit */ false);

  const int id = callbacks.add(static_callback);

  /** Id should be 1, 0 is loop_abort */
  Testing::reduce_and_check(world, id == 1);

  if (world.rank() == 0) {
    callbacks.call(static_callback, 5, 13);
    static_callback(5, 13);
    callbacks.abort_loop();
  } else {
    callbacks.loop();
  }

  Testing::reduce_and_check(world, static_callback_called);
}

/**
 * Check if adding and calling a dynamic callback
 * works as expected.
 */
BOOST_AUTO_TEST_CASE(add_dynamic_callback) {
  communicator world;
  MpiCallbacks callbacks(world, /* abort_on_exit */ false);

  bool called = false;
  auto cb = [&called](int a, int b) { called = (a == 5) && (b == 13); };

  const int id = callbacks.add(cb);

  /** Id should be 1, 0 is loop_abort */
  Testing::reduce_and_check(world, id == 1);

  if (world.rank() == 0) {
    callbacks.call(id, 5, 13);
    cb(5, 13);
    callbacks.abort_loop();
  } else {
    callbacks.loop();
  }

  Testing::reduce_and_check(world, called);
}

/**
 * Check wether remoing a dynamic callback
 * works.
 */
BOOST_AUTO_TEST_CASE(remove_dynamic_callback) {
  communicator world;
  MpiCallbacks callbacks(world, /* abort_on_exit */ false);

  bool called = false;
  auto cb = [&called](int a, int b) { called = (a == 5) && (b == 13); };

  const int dynamic_id = callbacks.add(cb);
  /** Id should be 1, 0 is loop_abort */
  Testing::reduce_and_check(world, dynamic_id == 1);

  const int static_id = callbacks.add(static_callback);
  /** Id should be 2 */
  Testing::reduce_and_check(world, static_id == 2);

  /** Remove dynamic callback */
  callbacks.remove(dynamic_id);

  if (world.rank() == 0) {
    /** Check that it is gone */
    BOOST_CHECK_THROW(callbacks.call(dynamic_id, 0, 0), std::out_of_range);
  }

  /** Calling the other callback should still work */
  static_callback_called = false;
  if (world.rank() == 0) {
    callbacks.call(static_callback, 5, 13);
    static_callback(5, 13);
    callbacks.abort_loop();
  } else {
    callbacks.loop();
  }

  Testing::reduce_and_check(world, static_callback_called);

  /** Re-adding should work. */
  const int new_dynamic_id = callbacks.add(cb);
  /** Id should be recycled id 1 */
  Testing::reduce_and_check(world, new_dynamic_id == 1);

  /** New callback should work */
  if (world.rank() == 0) {
    callbacks.call(new_dynamic_id, 5, 13);
    cb(5, 13);
    callbacks.abort_loop();
  } else {
    callbacks.loop();
  }

  Testing::reduce_and_check(world, called);
}

/* Check that the destructor calls abort */
BOOST_AUTO_TEST_CASE(destructor) {
  communicator world;

  {
    /* Will be detroyed on scope exit */
    MpiCallbacks callbacks(world);

    if (world.rank() == 0) {
      ;
    } else {
      callbacks.loop();
    }
  }

  /* This must be reachable by all nodes */
  Testing::reduce_and_check(world, true);
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
