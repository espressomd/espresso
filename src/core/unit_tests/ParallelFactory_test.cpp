/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
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

/** \file ParallelFactory_test.cpp Unit tests for the Utils::ParallelFactory class.
 *
*/

#include <string>

#include <boost/mpi.hpp>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ParallelFactory test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "utils/ParallelFactory.hpp"
#include "errorhandling.hpp"

namespace Testing {

void reduce_and_check(const boost::mpi::communicator &comm, bool local_value) {
  if(comm.rank() == 0) {
    bool total;
    boost::mpi::reduce(comm, local_value, total, std::logical_and<bool>(), 0);
    BOOST_CHECK(total);
  } else {
    boost::mpi::reduce(comm, local_value, std::logical_and<bool>(), 0);
  }
}

}

class TestClass {
 private:
  friend Utils::ParallelFactory<TestClass>;
  friend Utils::Factory<TestClass>;
  TestClass()
  {
    constructed = true;
  }
  ~TestClass()
  {
    destructed = true;
  }
 public:  
  static bool constructed, destructed;
};

bool TestClass::constructed{false};
bool TestClass::destructed{false};

using boost::mpi::communicator;
using Utils::ParallelFactory;
using std::string;

using Communication::mpiCallbacks;

BOOST_AUTO_TEST_CASE(make_instance) {
  const communicator &world = mpiCallbacks().comm();
  ParallelFactory<TestClass> factory;

  /* Register the test class with the underlying factory. */
  ParallelFactory<TestClass>::factory_type::register_new<TestClass>("TestClass");

  if(world.rank() == 0) {
    auto sp = factory.make("TestClass");
    mpiCallbacks().abort_loop();
    Testing::reduce_and_check(world, TestClass::constructed && sp->constructed);
  } else {
    mpiCallbacks().loop();
    Testing::reduce_and_check(world, TestClass::constructed);
  }

  if(world.rank() == 0) {
    BOOST_CHECK(check_runtime_errors() == 0);
  } else {
    check_runtime_errors();
  }
}

BOOST_AUTO_TEST_CASE(remove_instance) {
  const communicator &world = mpiCallbacks().comm();
  ParallelFactory<TestClass> factory;
  
  if(world.rank() == 0) {
    auto sp = factory.make("TestClass");

    /* Check that the reference count is checked. */
    {
      auto sp2 = sp;
      BOOST_CHECK_THROW(factory.remove(sp), std::runtime_error);
    }
    
    factory.remove(sp);

    /* sp should have been reset. */
    BOOST_CHECK(sp == nullptr);
    mpiCallbacks().abort_loop();
    Testing::reduce_and_check(world, TestClass::destructed);
  } else {
    mpiCallbacks().loop();
    Testing::reduce_and_check(world, TestClass::destructed);
  }

  if(world.rank() == 0) {
    BOOST_CHECK(check_runtime_errors() == 0);
  } else {
    check_runtime_errors();
  }
}

int main(int argc, char**argv) {  
  /* Init mpi */
  boost::mpi::environment mpi_env(argc, argv);
  boost::mpi::communicator world;
  
  /* Init callbacks */
  Communication::initialize_callbacks(world);
  /* Init error collector */
  ErrorHandling::init_error_handling(mpiCallbacks());

  /* Run tests */
  boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

