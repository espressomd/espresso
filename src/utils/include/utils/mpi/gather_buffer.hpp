/*
 * Copyright (C) 2017-2022 The ESPResSo project
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
 * The buffer is resized to the total size. On the @p root node,
 * the first @p n_elem elements of @p buffer are moved, if need
 * be. On the other nodes, @p buffer is not touched.
 *
 * This encapsulates a common combination of <tt>MPI_Gather()</tt>
 * and <tt>MPI_{Send,Recv}()</tt>.
 *
 * @param buffer On the head node: the target buffer that has the local
 *        part in the beginning. On worker nodes: the local buffer.
 * @param comm The MPI communicator.
 * @param root The rank where the data should be gathered.
 */
template <typename T, class Allocator>
void gather_buffer(std::vector<T, Allocator> &buffer,
                   boost::mpi::communicator const &comm, int root = 0) {
  auto const n_elem = static_cast<int>(buffer.size());

  if (comm.rank() == root) {
    static std::vector<int> sizes;
    static std::vector<int> displ;

    auto const tot_size =
        detail::size_and_offset<T>(sizes, displ, n_elem, comm, root);

    /* Resize the buffer */
    buffer.resize(static_cast<unsigned int>(tot_size));

    /* Move the original data to its new location */
    if (sizes[root] && displ[root]) {
      for (int i = sizes[root] - 1; i >= 0; --i) {
        buffer[i + displ[root]] = buffer[i];
      }
    }

    /* Gather data */
    gatherv(comm, buffer.data(), tot_size, buffer.data(), sizes.data(),
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
