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

#include "config/config.hpp"

#include "communication.hpp"

#include "cuda/init.hpp"
#include "errorhandling.hpp"
#include "fft/init.hpp"

#ifdef ESPRESSO_WALBERLA
#include <walberla_bridge/walberla_init.hpp>
#endif

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>
#include <omp.h>
#endif

#include <utils/Vector.hpp>
#include <utils/mpi/cart_comm.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <mpi.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <tuple>
#include <utility>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
namespace Communication {
struct KokkosHandle {
  KokkosHandle() { Kokkos::initialize(); }
  ~KokkosHandle() { Kokkos::finalize(); }
};
} // namespace Communication
#endif

boost::mpi::communicator comm_cart;
Communicator communicator{};
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
std::shared_ptr<Communication::KokkosHandle> kokkos_handle{};
#endif

namespace Communication {
static std::shared_ptr<MpiCallbacks> m_callbacks;

/* We use a singleton callback class for now. */
MpiCallbacks &mpiCallbacks() {
  assert(m_callbacks && "Mpi not initialized!");

  return *m_callbacks;
}

std::shared_ptr<MpiCallbacks> mpiCallbacksHandle() {
  assert(m_callbacks && "Mpi not initialized!");

  return m_callbacks;
}
} // namespace Communication

using Communication::mpiCallbacks;

int this_node = -1;

namespace Communication {
void init(std::shared_ptr<boost::mpi::environment> mpi_env) {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  const char *const env_omp_num_threads = std::getenv("OMP_NUM_THREADS");
  if (not env_omp_num_threads or std::strlen(env_omp_num_threads) == 0ul) {
    omp_set_num_threads(1);
  }
#endif

  communicator.full_initialization();

  Communication::m_callbacks =
      std::make_shared<Communication::MpiCallbacks>(comm_cart, mpi_env);

  ErrorHandling::init_error_handling(Communication::m_callbacks);

#ifdef ESPRESSO_WALBERLA
  walberla::mpi_init();
#endif

#ifdef ESPRESSO_CUDA
  cuda_on_program_start();
#endif

#ifdef ESPRESSO_FFTW
  fft_on_program_start();
#endif

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  kokkos_handle = std::make_shared<KokkosHandle>();
#endif
}

void deinit() {
  ErrorHandling::deinit_error_handling();
  Communication::m_callbacks.reset();

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  kokkos_handle.reset();
#endif
}
} // namespace Communication

Communicator::Communicator()
    : comm{::comm_cart}, node_grid{}, this_node{::this_node}, size{-1} {}

void Communicator::init_comm_cart() {
  auto constexpr reorder = false;
  comm = Utils::Mpi::cart_create(comm, node_grid, reorder);
  this_node = comm.rank();
  // check topology validity
  std::ignore = Utils::Mpi::cart_neighbors<3>(comm);
}

void Communicator::full_initialization() {
  assert(this_node == -1);
  assert(size == -1);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  node_grid = Utils::Mpi::dims_create<3>(size);
  init_comm_cart();
}

void Communicator::set_node_grid(Utils::Vector3i const &value) {
  node_grid = value;
  init_comm_cart();
}

Utils::Vector3i Communicator::calc_node_index() const {
  return Utils::Mpi::cart_coords<3>(comm, this_node);
}

std::shared_ptr<boost::mpi::environment> mpi_init(int argc, char **argv) {
  return std::make_shared<boost::mpi::environment>(argc, argv);
}

void mpi_loop() {
  if (this_node != 0)
    mpiCallbacks().loop();
}
