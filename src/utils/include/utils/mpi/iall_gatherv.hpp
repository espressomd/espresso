/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef UTILS_MPI_ALL_GATHERV_HPP
#define UTILS_MPI_ALL_GATHERV_HPP

#include "utils/Span.hpp"

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>

#include <algorithm>
#include <vector>

namespace Utils {
namespace Mpi {
namespace detail {

inline std::vector<int> displacements(Span<int const> sizes) {
  std::vector<int> displ(sizes.size());

  int offset = 0;
  for (int i = 0; i < displ.size(); i++) {
    displ[i] = offset;
    offset += sizes[i];
  }

  return displ;
}

template <typename T>
std::vector<boost::mpi::request>
iall_gatherv_impl(boost::mpi::communicator const &comm, T const *in_values,
                  int in_size, T *out_values, int const *sizes,
                  int const *displs) {
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

  return req;
}

} // namespace detail

template <typename T>
auto iall_gatherv(boost::mpi::communicator const &comm, T const *in_values,
                  int in_size, T *out_values, int const *sizes) {
  auto const displ =
      detail::displacements({sizes, static_cast<size_t>(comm.size())});

  return detail::iall_gatherv_impl(comm, in_values, in_size, out_values, sizes,
                                   displ.data());
}
} // namespace Mpi
} // namespace Utils
#endif
