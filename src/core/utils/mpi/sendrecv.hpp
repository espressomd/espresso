/*
  Copyright (C) 2018 The ESPResSo project

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

#ifndef UTILS_MPI_SENDRECV_HPP
#define UTILS_MPI_SENDRECV_HPP

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/mpi/nonblocking.hpp>

#include <array>

namespace Utils {
namespace Mpi {
namespace mpi = boost::mpi;

namespace detail {
template <typename T>
mpi::request isend(mpi::communicator const &comm, int dest, int tag,
                   const T &value, boost::mpl::false_) {
  /* boost.mpi before that version has a broken implementation
     of isend for non-mpi data types, so we have to use our own. */
#if BOOST_VERSION < 106800
  auto archive =
      boost::shared_ptr<mpi::packed_oarchive>(new mpi::packed_oarchive(comm));
  *archive << value;
  auto size = archive->size();
  mpi::request result;
  BOOST_MPI_CHECK_RESULT(MPI_Isend,
                         (&size, 1, mpi::get_mpi_datatype<std::size_t>(size),
                          dest, tag, comm, &result.m_requests[0]));
  BOOST_MPI_CHECK_RESULT(MPI_Isend,
                         (const_cast<void *>(archive->address()), size,
                          MPI_PACKED, dest, tag, comm, &result.m_requests[1]));

  result.m_data = archive;
  return result;
#else
  return comm.isend(dest, tag, value);
#endif
}

template <typename T>
mpi::request isend(mpi::communicator const &comm, int dest, int tag,
                   const T &value, boost::mpl::true_) {
  return comm.isend(dest, tag, value);
}

template <typename T>
std::array<mpi::request, 2> isendrecv_impl(mpi::communicator const &comm,
                                           int dest, int stag, const T &sval,
                                           int src, int rtag, T &rval) {
  return {{isend(comm, dest, stag, sval, mpi::is_mpi_datatype<T>()),
           comm.irecv(src, rtag, rval)}};
}

template <typename T>
mpi::status sendrecv_impl(mpi::communicator const &comm, int dest, int stag,
                          const T &sval, int src, int rtag, T &rval,
                          boost::mpl::true_) {
  mpi::status stat;
  BOOST_MPI_CHECK_RESULT(
      MPI_Sendrecv, (const_cast<T *>(&sval), 1, mpi::get_mpi_datatype<T>(sval),
                     dest, stag, &rval, 1, mpi::get_mpi_datatype<T>(rval), src,
                     rtag, comm, &stat.m_status));
  return stat;
}

template <typename T>
mpi::status sendrecv_impl(mpi::communicator const &comm, int dest, int stag,
                          const T &sval, int src, int rtag, T &rval,
                          boost::mpl::false_) {
  auto srrequests = isendrecv_impl(comm, dest, stag, sval, src, rtag, rval);
  mpi::status srstatuses[2];
  wait_all(srrequests.begin(), srrequests.end(), srstatuses);
  return srstatuses[1];
}
} // namespace detail

template <typename T>
mpi::status sendrecv(mpi::communicator const &comm, int dest, int stag,
                     const T &sval, int src, int rtag, T &rval) {
  return detail::sendrecv_impl(comm, dest, stag, sval, src, rtag, rval,
                               mpi::is_mpi_datatype<T>());
}

template <typename T>
std::array<mpi::request, 2> isendrecv(mpi::communicator const &comm, int dest,
                                      int stag, const T &sval, int src,
                                      int rtag, T &rval) {
  return detail::isendrecv_impl(comm, dest, stag, sval, src, rtag, rval);
}

} // namespace Mpi
} // namespace Utils

#endif // ESPRESSO_SENDRECV_HPP
