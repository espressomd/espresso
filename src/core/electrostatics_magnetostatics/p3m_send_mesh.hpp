//
// Created by florian on 17.11.19.
//

#ifndef ESPRESSO_P3M_SEND_MESH_HPP
#define ESPRESSO_P3M_SEND_MESH_HPP

#include "p3m-common.hpp"
#include <boost/mpi/communicator.hpp>

#include <vector>

/** Structure for send/recv meshes. */
struct p3m_send_mesh {
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

  void resize(const boost::mpi::communicator &comm,
              const p3m_local_mesh &local_mesh);
  void gather_grid(double *themesh, const boost::mpi::communicator &comm,
                   const int dim[3]);
  void spread_grid(double *themesh, const boost::mpi::communicator &comm,
                   const int dim[3]);
};

#endif // ESPRESSO_P3M_SEND_MESH_HPP
