/*
  Copyright (C) 2017 The ESPResSo project
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

#ifndef UTILS_MPI_DETAIL_SIZE_AND_OFFSET_HPP
#define UTILS_MPI_DETAIL_SIZE_AND_OFFSET_HPP

#include <algorithm>
#include <vector>
#include <numeric>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

namespace Utils {
namespace Mpi {
namespace detail {

template <typename T>
int size_and_offset(std::vector<int> &sizes, std::vector<int> &displ,
                    int n_elem, const boost::mpi::communicator &comm, int root = 0) {
  sizes.resize(comm.size());
  displ.resize(comm.size());

  /* Gather sizes */
  boost::mpi::gather(comm, n_elem, sizes, root);

  auto total_size = std::accumulate(sizes.begin(), sizes.end(), 0);

  int offset = 0;
  for (int i = 0; i < sizes.size(); i++) {
    displ[i] = offset;
    offset += sizes[i];
  }

  return total_size;
}

inline void size_and_offset(int n_elem, const boost::mpi::communicator &comm, int root = 0) {
  /* Send local size */
  boost::mpi::gather(comm, n_elem, root);
}
}
}
}

#endif
