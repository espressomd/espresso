/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

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
#if defined(OPEN_MPI)
#include <mpi-ext.h>
#endif

#include <cassert>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
struct KokkosHandle {
  KokkosHandle() { Kokkos::initialize(); }
  ~KokkosHandle() { Kokkos::finalize(); }
};
#endif

boost::mpi::communicator comm_cart;
Communicator communicator{};
std::unique_ptr<CommunicationEnvironment> communication_environment{};
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
std::shared_ptr<KokkosHandle> kokkos_handle{};
#endif
int this_node = -1;

[[maybe_unused]] static auto get_env_variable(char const *const name) {
  char const *const value = std::getenv(name);
  std::optional<std::string> result{std::nullopt};
  if (value) {
    result = std::string(value);
  }
  return result;
}

static auto make_default_mpi_env() {
  int argc = 0;
  char **argv = nullptr;
  return std::make_shared<boost::mpi::environment>(argc, argv);
}

CommunicationEnvironment::CommunicationEnvironment()
    : CommunicationEnvironment(make_default_mpi_env()) {}

CommunicationEnvironment::CommunicationEnvironment(
    std::shared_ptr<boost::mpi::environment> mpi_env)
    : m_mpi_env{std::move(mpi_env)} {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  auto const num_threads_env = get_env_variable("OMP_NUM_THREADS");
  if (not num_threads_env or num_threads_env->empty()) {
    omp_set_num_threads(1);
  }
#endif

  m_is_mpi_gpu_aware = false;

#if defined(OPEN_MPI)
#if defined(OMPI_HAVE_MPI_EXT_ROCM) && OMPI_HAVE_MPI_EXT_ROCM
  m_is_mpi_gpu_aware |= static_cast<bool>(MPIX_Query_rocm_support());
#endif
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
  m_is_mpi_gpu_aware |= static_cast<bool>(MPIX_Query_cuda_support());
#endif
#endif // defined(OPEN_MPI)

#if defined(MPICH)
  auto const mpich_gpu_env = get_env_variable("MPIR_CVAR_ENABLE_GPU");
  m_is_mpi_gpu_aware |= (mpich_gpu_env and *mpich_gpu_env == "1");
#endif // defined(MPICH)

#if defined(_CRAYC) or defined(__cray__)
  auto const cray_mpich_gpu_env = get_env_variable("MPICH_GPU_SUPPORT_ENABLED");
  m_is_mpi_gpu_aware |= (cray_mpich_gpu_env and *cray_mpich_gpu_env == "1");
#endif // defined(_CRAYC) or defined(__cray__)

#ifdef ESPRESSO_WALBERLA
  walberla::mpi_init();
#endif

  communicator.full_initialization();

  m_callbacks =
      std::make_shared<Communication::MpiCallbacks>(comm_cart, m_mpi_env);

  ErrorHandling::init_error_handling(comm_cart);

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

CommunicationEnvironment::~CommunicationEnvironment() {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  Kokkos::fence();
  kokkos_handle.reset();
#endif

#ifdef ESPRESSO_WALBERLA
  walberla::mpi_deinit();
#endif

  ErrorHandling::deinit_error_handling();
  m_callbacks.reset();
}

Communicator::Communicator()
    : comm{::comm_cart}, node_grid{}, this_node{::this_node}, size{-1},
      locked_for_checkpointing{false} {}

void Communicator::init_comm_cart() {
  auto constexpr reorder = false;
  comm = Utils::Mpi::cart_create(comm, node_grid, reorder);
  this_node = comm.rank();
  // check topology validity
  std::ignore = Utils::Mpi::cart_neighbors<3>(comm);
#ifdef ESPRESSO_WALBERLA
  walberla::mpi_reinit(node_grid.data());
#endif
}

void Communicator::full_initialization() {
  assert(this_node == -1);
  assert(size == -1);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  node_grid = Utils::Mpi::dims_create<3>(size);
  init_comm_cart();
}

void Communicator::set_node_grid(Utils::Vector3i const &value) {
  assert(not locked_for_checkpointing);
  node_grid = value;
  init_comm_cart();
}

Utils::Vector3i Communicator::calc_node_index() const {
  return Utils::Mpi::cart_coords<3>(comm, this_node);
}

void mpi_loop() {
  if (this_node != 0)
    Communication::mpiCallbacks().loop();
}
