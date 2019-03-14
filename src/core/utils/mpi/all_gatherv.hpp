/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef UTILS_MPI_ALL_GATHERV_HPP
#define UTILS_MPI_ALL_GATHERV_HPP

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <vector>

namespace Utils {
namespace Mpi {

namespace detail {
template <typename T>
void all_gatherv_impl(const boost::mpi::communicator &comm, const T *in_values,
                      int in_size, T *out_values, const int *sizes,
                      const int *displs, boost::mpl::true_) {
  MPI_Datatype type = boost::mpi::get_mpi_datatype<T>();

  /* in-place ? */
  if (in_values == out_values) {
    BOOST_MPI_CHECK_RESULT(MPI_Allgatherv,
                           (MPI_IN_PLACE, 0, type, out_values,
                            const_cast<int *>(sizes), const_cast<int *>(displs),
                            type, comm));
  } else {
    BOOST_MPI_CHECK_RESULT(MPI_Allgatherv,
                           (const_cast<T *>(in_values), in_size, type,
                            out_values, const_cast<int *>(sizes),
                            const_cast<int *>(displs), type, comm));
  }
}

template <typename T>
void all_gatherv_impl(const boost::mpi::communicator &comm, const T *in_values,
                      int in_size, T *out_values, const int *sizes,
                      const int *displs, boost::mpl::false_) {
  auto const n_nodes = comm.size();
  auto const rank = comm.rank();

  /* not in-place */
  if (in_values != out_values) {
    std::copy_n(in_values, in_size, out_values + displs[rank]);
  }

  std::vector<boost::mpi::request> req;
  for (int i = 0; i < n_nodes; i++) {
    if (i != rank) {
      req.emplace_back(comm.isend(i, 42, out_values + displs[rank], in_size));
      req.emplace_back(comm.irecv(i, 42, out_values + displs[i], sizes[i]));
    }
  }

  boost::mpi::wait_all(req.begin(), req.end());
}
} // namespace detail

template <typename T>
void all_gatherv(const boost::mpi::communicator &comm, const T *in_values,
                 int in_size, T *out_values, const int *sizes,
                 const int *displs) {
  detail::all_gatherv_impl(comm, in_values, in_size, out_values, sizes, displs,
                           boost::mpi::is_mpi_datatype<T>());
}

template <typename T>
void all_gatherv(const boost::mpi::communicator &comm, const T *in_values,
                 int in_size, T *out_values, const int *sizes) {
  std::vector<int> displ(comm.size());

  int offset = 0;
  for (unsigned i = 0; i < displ.size(); i++) {
    displ[i] = offset;
    offset += sizes[i];
  }

  detail::all_gatherv_impl(comm, in_values, in_size, out_values, sizes,
                           displ.data(), boost::mpi::is_mpi_datatype<T>());
}
} // namespace Mpi
} // namespace Utils
#endif
