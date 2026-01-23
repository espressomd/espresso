/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#include <walberla_bridge/utils/ResourceManager.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <core/mpi/Environment.h>
#include <core/mpi/MPIManager.h>

#include <cassert>
#include <memory>

/** @brief waLBerla MPI communicator. */
static std::shared_ptr<walberla::mpi::MPIManager> walberla_mpi_comm;
/** @brief waLBerla MPI environment. */
static std::shared_ptr<walberla::mpi::Environment> walberla_mpi_env;
/**
 * @brief waLBerla MPI Cartesian communicator observer.
 * It is purposefully decoupled from @ref walberla_mpi_comm.
 * We are not tracking the lifetime of @ref walberla_mpi_comm itself,
 * but the lifetime of its ownership by the waLBerla MPI singleton.
 */
static std::shared_ptr<int> walberla_mpi_cart_comm_observer;

namespace walberla {

void mpi_init() {
  assert(::walberla_mpi_env == nullptr);
  assert(::walberla_mpi_comm == nullptr);
  assert(::walberla_mpi_cart_comm_observer == nullptr);
  int argc = 0;
  char **argv = nullptr;
  ::walberla_mpi_env = std::make_shared<walberla::mpi::Environment>(argc, argv);
  ::walberla_mpi_comm = walberla::mpi::MPIManager::instance();
  ::walberla_mpi_cart_comm_observer = std::make_shared<int>(0);
}

void mpi_reinit(int const *const cart_topol) {
  assert(::walberla_mpi_env.use_count() >= 1);
  assert(::walberla_mpi_comm.use_count() >= 1);
  assert(::walberla_mpi_cart_comm_observer.use_count() == 1);
  ::walberla_mpi_comm->resetMPI();
  ::walberla_mpi_comm->createCartesianComm(cart_topol[0], cart_topol[1],
                                           cart_topol[2], true, true, true);
  ::walberla_mpi_cart_comm_observer = std::make_shared<int>(0);
}

void mpi_deinit() {
  assert(::walberla_mpi_env.use_count() >= 1);
  assert(::walberla_mpi_comm.use_count() >= 1);
  assert(::walberla_mpi_cart_comm_observer.use_count() == 1);
  ::walberla_mpi_env.reset();
  ::walberla_mpi_comm.reset();
  ::walberla_mpi_cart_comm_observer.reset();
}

std::unique_ptr<ResourceManager> get_vtk_dependent_resources() {
  auto vtk_dependencies = std::make_unique<ResourceManager>();
  // waLBerla MPI communicator (singleton)
  vtk_dependencies->acquire_lock(::walberla_mpi_comm);
  // waLBerla MPI environment (destructor depends on the MPI communicator)
  vtk_dependencies->acquire_lock(::walberla_mpi_env);
  return vtk_dependencies;
}

ResourceObserver get_mpi_cart_comm_observer() {
  assert(walberla_mpi_cart_comm_observer.use_count() <= 1);
  return walberla_mpi_cart_comm_observer;
}

} // namespace walberla
