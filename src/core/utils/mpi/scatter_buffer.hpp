/*
  Copyright (C) 2017-2018 The ESPResSo project
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

#ifndef UTILS_MPI_SCATTER_BUFFER_HPP
#define UTILS_MPI_SCATTER_BUFFER_HPP

#include "detail/size_and_offset.hpp"

#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <type_traits>
#include <vector>

namespace Utils {
namespace Mpi {

/**
 * @brief Scatter buffer with different size on each node.
 *
 * Scatter a buffer to the nodes, where every node gets
 * a different chunk of the buffer, controlled by the slave.
 *
 * This is a collective call.
 */
template <typename T>
void scatter_buffer(T *buffer, int n_elem, boost::mpi::communicator comm,
                    int root = 0) {
  static_assert(std::is_pod<T>::value, "");
  if (comm.rank() == root) {
    static std::vector<int> sizes;
    static std::vector<int> displ;

    detail::size_and_offset<T>(sizes, displ, n_elem, comm, root);

    for (int i = 0; i < comm.size(); i++) {
      sizes[i] *= sizeof(T);
      displ[i] *= sizeof(T);
    }

    /* Send data */
    MPI_Scatterv(buffer, sizes.data(), displ.data(), MPI_BYTE, MPI_IN_PLACE, 0,
                 MPI_BYTE, root, comm);
  } else {
    detail::size_and_offset(n_elem, comm, root);
    /* Recv data */
    MPI_Scatterv(nullptr, nullptr, nullptr, MPI_BYTE, buffer,
                 n_elem * sizeof(T), MPI_BYTE, root, comm);
  }
}
} // namespace Mpi
} // namespace Utils

#endif
