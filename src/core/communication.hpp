/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/** \file
 *  This file contains the asynchronous MPI communication.
 *
 *  It is the header file for communication.cpp.
 *
 *  The asynchronous MPI communication is used during the script
 *  evaluation. Except for the head node that interprets the interface
 *  script, all other nodes wait in @ref mpi_loop() for the head node to
 *  issue an action using @c MpiCallbacks::call(). @ref mpi_loop() immediately
 *  executes an @c MPI_Bcast and therefore waits for the head node to
 *  broadcast a command, which is done by @c MpiCallbacks::call(). The request
 *  consists of a callback function with an arbitrary number of arguments.
 *
 *  To add new actions (e.g. to implement new interface functionality), do the
 *  following:
 *  - write the @c mpi_* function that is executed on the head node
 *  - write the @c mpi_*_local function that is executed on worker nodes
 *  - register the local function with one of the @c REGISTER_CALLBACK macros
 *
 *  After this, your procedure is free to do anything. However, it has
 *  to be in (MPI) sync with what your new @c mpi_*_local does. This
 *  procedure is called immediately after the broadcast.
 */

#include "MpiCallbacks.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <memory>
#include <utility>

/** The number of this node. */
extern int this_node;
/** The communicator */
extern boost::mpi::communicator comm_cart;
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
struct KokkosHandle;
extern std::shared_ptr<KokkosHandle> kokkos_handle;
#endif

class CommunicationEnvironment {
  std::shared_ptr<boost::mpi::environment> m_mpi_env;
  std::shared_ptr<Communication::MpiCallbacks> m_callbacks;
  bool m_is_mpi_gpu_aware;

public:
  CommunicationEnvironment();
  explicit CommunicationEnvironment(
      std::shared_ptr<boost::mpi::environment> mpi_env);
  ~CommunicationEnvironment();

  auto &mpiCallbacks() const { return *m_callbacks; }
  auto mpiCallbacksHandle() { return m_callbacks; }
  auto get_mpi_env() const { return m_mpi_env; }
  auto is_mpi_gpu_aware() const { return m_is_mpi_gpu_aware; }
};

struct Communicator {
  boost::mpi::communicator &comm;
  Utils::Vector3i node_grid;
  /** @brief The MPI rank. */
  int &this_node;
  /** @brief The MPI world size. */
  int size;

  Communicator();
  void init_comm_cart();
  void full_initialization();
  /** @brief Calculate the node index in the Cartesian topology. */
  Utils::Vector3i calc_node_index() const;
  /** @brief Set new Cartesian topology. */
  void set_node_grid(Utils::Vector3i const &value);
};

extern Communicator communicator;
extern std::unique_ptr<CommunicationEnvironment> communication_environment;

namespace Communication {
/**
 * @brief Returns a reference to the global callback class instance.
 */
inline MpiCallbacks &mpiCallbacks() {
  return ::communication_environment->mpiCallbacks();
}
} // namespace Communication

/** Process requests from head node. Worker nodes main loop. */
void mpi_loop();
