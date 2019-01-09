// Copyright 2005 Douglas Gregor.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Message Passing Interface 1.1 -- Section 3. MPI Point-to-point

/* There is the potential for optimization here. We could keep around
   a "small message" buffer of size N that we just receive into by
   default. If the message is N - sizeof(int) bytes or smaller, it can
   just be sent with that buffer. If it's larger, we send the first N
   - sizeof(int) bytes in the first packet followed by another
   packet. The size of the second packet will be stored in an integer
   at the end of the first packet.

   We will introduce this optimization later, when we have more
   performance test cases and have met our functionality goals. */

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/mpi/detail/point_to_point.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/mpi/request.hpp>

template <>
boost::mpi::request
boost::mpi::communicator::isend<boost::mpi::packed_oarchive>(
    int dest, int tag, boost::mpi::packed_oarchive const &ar) const {
  std::size_t const& size = ar.size();
  request req;
  BOOST_MPI_CHECK_RESULT(MPI_Isend,
                         (&const_cast<std::size_t &>(size), 1,
                                 get_mpi_datatype(size),
                                 dest, tag, *this, &req.m_requests[0]));
  BOOST_MPI_CHECK_RESULT(MPI_Isend,
                         (const_cast<void *>(ar.address()), size,
                                 MPI_PACKED,
                                 dest, tag, *this, &req.m_requests[1]));

  return req;
}
