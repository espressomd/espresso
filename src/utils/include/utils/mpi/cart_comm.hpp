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
#ifndef UTILS_MPI_CART_COMM_HPP
#define UTILS_MPI_CART_COMM_HPP

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/exception.hpp>

#include <mpi.h>

namespace Utils {
namespace Mpi {

/**
 * @brief Wrapper around MPI_Dims_create.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dim> Vector<int, dim> dims_create(int nodes) {
  Vector<int, dim> dims{};
  BOOST_MPI_CHECK_RESULT(MPI_Dims_create,
                         (nodes, static_cast<int>(dim), dims.data()))

  return dims;
}

/**
 * @brief Wrapper around MPI_Cart_create.
 *
 * @tparam dim Number of dimensions
 */
template <size_t dim>
boost::mpi::communicator cart_create(
    boost::mpi::communicator const &comm, Vector<int, dim> const &dims,
    bool reorder = true,
    Vector<int, dim> const &periodicity = Vector<int, dim>::broadcast(1)) {
  MPI_Comm temp_comm;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_create,
                         (comm, dim, dims.data(), periodicity.data(),
                          static_cast<int>(reorder), &temp_comm))

  return boost::mpi::communicator(temp_comm, boost::mpi::comm_take_ownership);
}

/**
 * @brief Wrapper around MPI_Cart_coords.
 *
 * @tparam dims Number of dimensions
 */
template <size_t dims>
Vector3i cart_coords(boost::mpi::communicator const &comm, int rank) {
  Vector3i pos;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_coords, (comm, rank, dims, pos.data()))
  return pos;
}

/**
 * @brief Wrapper around MPI_Cart_rank.
 *
 * @tparam dims Number of dimensions
 */
template <size_t dims>
int cart_rank(boost::mpi::communicator const &comm,
              const Vector<int, dims> &pos) {
  int rank;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_rank, (comm, pos.data(), &rank))
  return rank;
}

/**
 * @brief Wrapper around MPI_Cart_shift.
 *
 * @return pair of source and destination rank.
 */
inline std::pair<int, int> cart_shift(boost::mpi::communicator const &comm,
                                      int direction, int displacement) {
  int src = -1, dst = -1;
  BOOST_MPI_CHECK_RESULT(MPI_Cart_shift,
                         (comm, direction, displacement, &src, &dst))

  return {src, dst};
}
} // namespace Mpi
} // namespace Utils

#endif
