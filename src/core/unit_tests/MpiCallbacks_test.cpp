#define BOOST_TEST_MODULE NumeratedContainerTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>

#include "MpiCallbacks.hpp"

boost::mpi::environment mpi_env;

template<int id>
void callback(int a, int b) {
  bool res = false;

  switch(id) {
    case 0:
      /** Check parameters */
      res = (a == 1) and (b == 2);
      
      break;
    case 1:
      /** Terminate */
      std::exit(0);
      
      break;
  }
  
  boost::mpi::reduce(MpiCallbacks::mpi_comm, res, res, std::logical_and<bool>(), 0);
}

BOOST_AUTO_TEST_CASE(test_MpiCallbacks) {
  BOOST_CHECK(MpiCallbacks::add(callback<0>) == 0);

  if(MpiCallbacks::mpi_comm.rank() == 0) {
      BOOST_CHECK_THROW(MpiCallbacks::call(1, 0, 0), std::out_of_range);
      BOOST_CHECK_THROW(MpiCallbacks::call(nullptr, 0, 0), std::out_of_range);
  }

  BOOST_CHECK(MpiCallbacks::add(callback<1>) == 1);
  
  if(MpiCallbacks::mpi_comm.rank() == 0) {
    MpiCallbacks::call(0, 1, 2);

    bool res = true, glob_res = true;
    boost::mpi::reduce(MpiCallbacks::mpi_comm, res, glob_res, std::logical_and<bool>(), 0);
    BOOST_CHECK(glob_res);

    MpiCallbacks::call(callback<0>, 1, 2);

    res = true;
    glob_res = true;
    boost::mpi::reduce(MpiCallbacks::mpi_comm, res, glob_res, std::logical_and<bool>(), 0);
    BOOST_CHECK(glob_res);
        
    /** Terminate slaves */
    MpiCallbacks::call(1, 0, 0);
  } else {
    MpiCallbacks::loop();
  }  
}
