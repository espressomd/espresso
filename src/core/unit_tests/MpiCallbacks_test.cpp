/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

/* Unit tests for the MpiCallbacks class. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE MpiCallbacks test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../MpiCallbacks.hpp"

#include <boost/mpi.hpp>
#include <boost/optional.hpp>

#include <string>

static bool called = false;

BOOST_AUTO_TEST_CASE(invoke_test) {
  using Communication::detail::invoke;

  auto f = [](int i, unsigned j) { return i + j; };

  boost::mpi::communicator world;
  boost::mpi::packed_oarchive::buffer_type buff;

  auto const i = 123;
  auto const j = 456u;
  boost::mpi::packed_oarchive(world, buff) << i << j;

  boost::mpi::packed_iarchive ia(world, buff);

  BOOST_CHECK_EQUAL(f(i, j), (invoke<decltype(f), int, unsigned>(f, ia)));
}

/*
 * Test that the implementation of callback_model_t
 * correctly deserialize the parameters and call
 * the callback with them.
 */
BOOST_AUTO_TEST_CASE(callback_model_t) {
  using namespace Communication;
  boost::mpi::communicator world;

  boost::mpi::packed_oarchive::buffer_type buff;
  boost::mpi::packed_oarchive oa(world, buff);

  {
    int i = 537;
    double d = 3.4;
    oa << i << d;
  }

  /* function pointer variant */
  {
    called = false;
    void (*fp)(int, double) = [](int i, double d) {
      BOOST_CHECK_EQUAL(537, i);
      BOOST_CHECK_EQUAL(3.4, d);

      called = true;
    };

    auto cb = detail::make_model(fp);

    boost::mpi::packed_iarchive ia(world, buff);
    cb->operator()(world, ia);

    BOOST_CHECK(called);
  }

  /* Lambda */
  {
    called = false;
    auto cb = detail::make_model([state = 19](int i, double d) {
      BOOST_CHECK_EQUAL(19, state);
      BOOST_CHECK_EQUAL(537, i);
      BOOST_CHECK_EQUAL(3.4, d);

      called = true;
    });

    boost::mpi::packed_iarchive ia(world, buff);
    cb->operator()(world, ia);

    BOOST_CHECK(called);
  }
}

BOOST_AUTO_TEST_CASE(adding_function_ptr_cb) {
  boost::mpi::communicator world;
  Communication::MpiCallbacks cb(world);

  void (*fp)(int, const std::string &) = [](int i, const std::string &s) {
    BOOST_CHECK_EQUAL(537, i);
    BOOST_CHECK_EQUAL("adding_function_ptr_cb", s);

    called = true;
  };

  cb.add(fp);

  called = false;

  if (0 == world.rank()) {
    cb.call(fp, 537, std::string("adding_function_ptr_cb"));
  } else {
    cb.loop();
    BOOST_CHECK(called);
  }
}

BOOST_AUTO_TEST_CASE(RegisterCallback) {
  void (*fp)(int, const std::string &) = [](int i, const std::string &s) {
    BOOST_CHECK_EQUAL(537, i);
    BOOST_CHECK_EQUAL("2nd", s);

    called = true;
  };

  Communication::RegisterCallback{fp};

  boost::mpi::communicator world;
  Communication::MpiCallbacks cb(world);

  called = false;

  if (0 == world.rank()) {
    cb.call(fp, 537, std::string("2nd"));
  } else {
    cb.loop();
    BOOST_CHECK(called);
  }
}

BOOST_AUTO_TEST_CASE(CallbackHandle) {
  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  bool m_called = false;
  Communication::CallbackHandle<std::string> cb(
      &cbs, [&m_called](std::string s) {
        BOOST_CHECK_EQUAL("CallbackHandle", s);

        m_called = true;
      });

  if (0 == world.rank()) {
    cb(std::string("CallbackHandle"));
  } else {
    cbs.loop();
    BOOST_CHECK(called);
  }
}

BOOST_AUTO_TEST_CASE(reduce_callback) {
  auto cb = []() -> int { return boost::mpi::communicator().rank(); };
  Communication::MpiCallbacks::add_static(Communication::Result::Reduction{},
                                          static_cast<int (*)()>(cb),
                                          std::plus<int>());

  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  if (0 == world.rank()) {
    auto const ret = cbs.call(Communication::Result::reduction,
                              std::plus<int>(), static_cast<int (*)()>(cb));
    auto const n = world.size();
    BOOST_CHECK_EQUAL(ret, (n * (n - 1)) / 2);
  } else {
    cbs.loop();
  }
}

BOOST_AUTO_TEST_CASE(ignore_callback) {
  called = false;
  auto cb = []() -> int { return called = true, -1; };

  auto const fp = static_cast<int (*)()>(cb);

  Communication::MpiCallbacks::add_static(Communication::Result::ignore, fp);

  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  if (0 == world.rank()) {
    cbs.call(Communication::Result::ignore, fp);
  } else {
    cbs.loop();
  }

  BOOST_CHECK(called);
}

BOOST_AUTO_TEST_CASE(one_rank_callback) {
  auto cb = []() -> boost::optional<int> {
    boost::mpi::communicator world;
    if (world.rank() == (world.size() - 1)) {
      return world.rank();
    }

    return {};
  };

  auto const fp = static_cast<boost::optional<int> (*)()>(cb);

  Communication::MpiCallbacks::add_static(Communication::Result::one_rank, fp);

  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  if (0 == world.rank()) {
    BOOST_CHECK_EQUAL(cbs.call(Communication::Result::one_rank, fp),
                      world.size() - 1);
  } else {
    cbs.loop();
  }
}

BOOST_AUTO_TEST_CASE(master_rank_callback) {
  auto cb = []() -> int {
    boost::mpi::communicator world;
    if (world.rank() == 0) {
      return world.size();
    }

    return -1;
  };

  auto const fp = static_cast<int (*)()>(cb);

  Communication::MpiCallbacks::add_static(Communication::Result::master_rank,
                                          fp);

  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  if (0 == world.rank()) {
    BOOST_CHECK_EQUAL(cbs.call(Communication::Result::master_rank, fp),
                      world.size());
  } else {
    cbs.loop();
  }
}

BOOST_AUTO_TEST_CASE(call_all) {
  called = false;
  auto cb = []() { called = true; };

  auto const fp = static_cast<void (*)()>(cb);

  Communication::MpiCallbacks::add_static(fp);

  boost::mpi::communicator world;
  Communication::MpiCallbacks cbs(world);

  if (0 == world.rank()) {
    cbs.call_all(fp);
  } else {
    cbs.loop();
  }

  BOOST_CHECK(called);
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
