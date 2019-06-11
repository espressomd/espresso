/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#ifndef UTILS_MPI_ALL_COMPARED_HPP
#define UTILS_MPI_ALL_COMPARED_HPP

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

namespace Utils {
namespace Mpi {

/**
 * @brief Compare values on all nodes.
 *
 * @param comm The communicator to operate on.
 * @param value the value to compare.
 *
 * @tparam T the type of value.
 *
 * Returns true iff all ranks in comm called
 * the function with the same value.
 * T has to be serializeable.
 */
template <typename T>
bool all_compare(boost::mpi::communicator const &comm, T const &value) {
  T root_value{};

  /* Arbitrary pick the value on rank 0, and broadcast it */
  if (comm.rank() == 0) {
    root_value = value;
  }

  boost::mpi::broadcast(comm, root_value, 0);

  bool is_same = false;
  boost::mpi::all_reduce(comm, value == root_value, is_same,
                         std::logical_and<bool>());

  return is_same;
}
} // namespace Mpi
} // namespace Utils

#endif
