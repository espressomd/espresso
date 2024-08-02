/*
 * Copyright (C) 2010-2024 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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

#pragma once

#include "config/config.hpp"

#if defined(P3M) or defined(DP3M)

#include "p3m/common.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <span>
#include <vector>

/** @brief P3M halo communicator. */
template <typename FloatType> class p3m_send_mesh {
  enum Requests {
    REQ_P3M_INIT = 200,
    REQ_P3M_GATHER = 201,
    REQ_P3M_SPREAD = 202
  };
  /** dimension of sub meshes to send. */
  int s_dim[6][3];
  /** lower left corners of sub meshes to send. */
  int s_ld[6][3];
  /** upper right corners of sub meshes to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimension of sub meshes to recv. */
  int r_dim[6][3];
  /** lower left corners of sub meshes to recv. */
  int r_ld[6][3];
  /** upper right corners of sub meshes to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;

  /** vector to store grid points to send. */
  std::vector<FloatType> send_grid;
  /** vector to store grid points to recv */
  std::vector<FloatType> recv_grid;

public:
  void resize(boost::mpi::communicator const &comm,
              P3MLocalMesh const &local_mesh);
  void gather_grid(boost::mpi::communicator const &comm,
                   std::span<FloatType *> meshes, Utils::Vector3i const &dim);
  void gather_grid(boost::mpi::communicator const &comm, FloatType *mesh,
                   Utils::Vector3i const &dim) {
    gather_grid(comm, std::span(&mesh, 1u), dim);
  }
  void spread_grid(boost::mpi::communicator const &comm,
                   std::span<FloatType *> meshes, Utils::Vector3i const &dim);
  void spread_grid(boost::mpi::communicator const &comm, FloatType *mesh,
                   Utils::Vector3i const &dim) {
    spread_grid(comm, std::span(&mesh, 1u), dim);
  }
};

#endif // defined(P3M) or defined(DP3M)
