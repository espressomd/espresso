/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_MPI_GATHERV_HPP
#define UTILS_MPI_GATHERV_HPP

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <vector>

namespace Utils {
namespace Mpi {

namespace detail {
template <typename T>
void gatherv_impl(const boost::mpi::communicator &comm, const T *in_values,
                  int in_size, T *out_values, const int *sizes,
                  const int *displs, int root, boost::mpl::true_) {
  if (in_values == nullptr)
    return;

  MPI_Datatype type = boost::mpi::get_mpi_datatype<T>(*in_values);

  /* in-place ? */
  if ((in_values == out_values) && (comm.rank() == root)) {
    BOOST_MPI_CHECK_RESULT(MPI_Gatherv,
                           (MPI_IN_PLACE, 0, type, out_values,
                            const_cast<int *>(sizes), const_cast<int *>(displs),
                            type, root, comm));
  } else {
    BOOST_MPI_CHECK_RESULT(MPI_Gatherv,
                           (const_cast<T *>(in_values), in_size, type,
                            out_values, const_cast<int *>(sizes),
                            const_cast<int *>(displs), type, root, comm));
  }
}

template <typename T>
void gatherv_impl(const boost::mpi::communicator &comm, const T *in_values,
                  int in_size, T *out_values, const int *sizes,
                  const int *displs, int root, boost::mpl::false_) {
  if (comm.rank() == root) {
    auto const n_nodes = comm.size();

    /* not in-place */
    if (in_values && out_values && in_values != out_values) {
      std::copy_n(in_values, in_size, out_values + displs[root]);
    }

    std::vector<boost::mpi::request> req;
    for (int i = 0; i < n_nodes; i++) {
      if (i == root)
        continue;

      // NOLINTNEXTLINE(clang-analyzer-core.NullDereference)
      req.emplace_back(comm.irecv(i, 42, out_values + displs[i], sizes[i]));
    }

    boost::mpi::wait_all(req.begin(), req.end());
  } else {
    comm.send(root, 42, in_values, in_size);
  }
}
} // namespace detail

template <typename T>
void gatherv(const boost::mpi::communicator &comm, const T *in_values,
             int in_size, T *out_values, const int *sizes, const int *displs,
             int root) {
  detail::gatherv_impl(comm, in_values, in_size, out_values, sizes, displs,
                       root, boost::mpi::is_mpi_datatype<T>());
}

template <typename T>
void gatherv(const boost::mpi::communicator &comm, const T *in_values,
             int in_size, T *out_values, const int *sizes, int root) {
  if (comm.rank() == root) {
    std::vector<int> displ(comm.size());

    int offset = 0;
    for (unsigned i = 0; i < displ.size(); i++) {
      displ[i] = offset;
      offset += sizes[i];
    }

    detail::gatherv_impl(comm, in_values, in_size, out_values, sizes,
                         displ.data(), root, boost::mpi::is_mpi_datatype<T>());

  } else {
    detail::gatherv_impl(comm, in_values, in_size, out_values, nullptr, nullptr,
                         root, boost::mpi::is_mpi_datatype<T>());
  }
}

template <typename T>
void gatherv(const boost::mpi::communicator &comm, const T *in_values,
             int in_size, int root) {
  assert(comm.rank() != root &&
         "This overload can not be called on the root rank.");
  gatherv(comm, in_values, in_size, static_cast<T *>(nullptr), nullptr, nullptr,
          root);
}

} // namespace Mpi
} // namespace Utils
#endif
