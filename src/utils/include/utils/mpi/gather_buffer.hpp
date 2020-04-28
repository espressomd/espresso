/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#ifndef UTILS_MPI_GATHER_BUFFER_HPP
#define UTILS_MPI_GATHER_BUFFER_HPP

#include "detail/size_and_offset.hpp"
#include "gatherv.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <type_traits>
#include <vector>

namespace Utils {
namespace Mpi {

/**
 * @brief Gather buffer with different size on each node.
 *
 * Gathers buffers with different lengths from all nodes to root.
 * The buffer is assumed to be large enough to hold the data from
 * all the nodes and is owned by the caller. On the root node no
 * data is copied, and the first n_elem elements of buffer are not
 * touched. This combines a common combination of MPI_Gather and
 * MPI_{Send,Recv}.
 *
 * @param buffer On the master the target buffer that has to be
          large enough to hold all elements and has the local
          part in the beginning. On the slaves the local buffer.
 * @param n_elem The number of elements in the local buffer.
 * @param comm The MPI communicator.
 * @param root The rank where the data should be gathered.
 * @return On rank root, the total number of elements in the buffer,
 *         on the other ranks 0.
 */
template <typename T>
int gather_buffer(T *buffer, int n_elem, boost::mpi::communicator comm,
                  int root = 0) {
  if (comm.rank() == root) {
    static std::vector<int> sizes;
    static std::vector<int> displ;

    auto const total_size =
        detail::size_and_offset<T>(sizes, displ, n_elem, comm, root);

    /* Gather data */
    gatherv(comm, buffer, 0, buffer, sizes.data(), displ.data(), root);

    return total_size;
  }
  detail::size_and_offset(n_elem, comm, root);
  /* Send data */
  gatherv(comm, buffer, n_elem, static_cast<T *>(nullptr), nullptr, nullptr,
          root);

  return 0;
}

/**
 * @brief Gather buffer with different size on each node.
 *
 * Gathers buffers with different lengths from all nodes to root.
 * The buffer is resized to the total size. On the root node no
 * data is copied, and the first n_elem elements of buffer are not
 * touched. On the slaves, the buffer is not touched.
 *
 * @param buffer On the master the target buffer that has the local
          part in the beginning. On the slaves the local buffer.
 * @param comm The MPI communicator.
 * @param root The rank where the data should be gathered.
 */
template <typename T, class Allocator>
void gather_buffer(std::vector<T, Allocator> &buffer,
                   boost::mpi::communicator comm, int root = 0) {
  auto const n_elem = buffer.size();

  if (comm.rank() == root) {
    static std::vector<int> sizes;
    static std::vector<int> displ;

    auto const tot_size =
        detail::size_and_offset<T>(sizes, displ, n_elem, comm, root);

    /* Resize the buffer */
    buffer.resize(tot_size);

    /* Gather data */
    gatherv(comm, buffer.data(), buffer.size(), buffer.data(), sizes.data(),
            displ.data(), root);
  } else {
    /* Send local size */
    detail::size_and_offset(n_elem, comm, root);
    /* Send data */
    gatherv(comm, buffer.data(), n_elem, static_cast<T *>(nullptr), nullptr,
            nullptr, root);
  }
}
} // namespace Mpi
} // namespace Utils

#endif
