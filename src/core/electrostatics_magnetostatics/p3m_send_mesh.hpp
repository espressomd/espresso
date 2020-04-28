/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ESPRESSO_P3M_SEND_MESH_HPP
#define ESPRESSO_P3M_SEND_MESH_HPP

#include "p3m-common.hpp"

#if defined(P3M) || defined(DP3M)

#include <utils/Span.hpp>

#include <boost/mpi/communicator.hpp>

#include <vector>

/** Structure for send/recv meshes. */
class p3m_send_mesh {
  enum Requests {
    REQ_P3M_INIT = 200,
    REQ_P3M_GATHER = 201,
    REQ_P3M_SPREAD = 202
  };
  /** dimension of sub meshes to send. */
  int s_dim[6][3];
  /** left down corners of sub meshes to send. */
  int s_ld[6][3];
  /** up right corners of sub meshes to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimension of sub meshes to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshes to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshes to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;

  /** vector to store grid points to send. */
  std::vector<double> send_grid;
  /** vector to store grid points to recv */
  std::vector<double> recv_grid;

public:
  void resize(const boost::mpi::communicator &comm,
              const p3m_local_mesh &local_mesh);
  void gather_grid(Utils::Span<double *> meshes,
                   const boost::mpi::communicator &comm,
                   const Utils::Vector3i &dim);
  void gather_grid(double *mesh, const boost::mpi::communicator &comm,
                   const Utils::Vector3i &dim) {
    gather_grid(Utils::make_span(&mesh, 1), comm, dim);
  }
  void spread_grid(Utils::Span<double *> meshes,
                   const boost::mpi::communicator &comm,
                   const Utils::Vector3i &dim);
  void spread_grid(double *mesh, const boost::mpi::communicator &comm,
                   const Utils::Vector3i &dim) {
    spread_grid(Utils::make_span(&mesh, 1), comm, dim);
  }
};
#endif
#endif // ESPRESSO_P3M_SEND_MESH_HPP
