/*
 * Copyright (C) 2026 The ESPResSo project
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

#ifdef ESPRESSO_P3M

#include "p3m/P3MFFTBackendFactory.hpp"

#ifdef ESPRESSO_KOKKOS_FFT
#include "p3m/P3MFFTKokkos.hpp"
#endif

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <memory>

template <typename FloatType>
std::shared_ptr<P3MFFTBackend<FloatType, P3MFFTKokkosConfig>>
make_p3m_kokkos_fft_backend(
    [[maybe_unused]] boost::mpi::communicator const &comm,
    [[maybe_unused]] Utils::Vector3i const &global_mesh,
    [[maybe_unused]] Utils::Vector3i const &rs_local_ld_index,
    [[maybe_unused]] Utils::Vector3i const &rs_local_ur_index,
    [[maybe_unused]] Utils::Vector3i const &node_grid) {
#ifdef ESPRESSO_KOKKOS_FFT
  // kokkos-fft is a local (non-MPI) transform: use it only on a single rank.
  if (comm.size() == 1) {
    return std::make_shared<P3MFFTKokkos<FloatType, P3MFFTKokkosConfig>>(
        comm, global_mesh, rs_local_ld_index, rs_local_ur_index, node_grid);
  }
#endif
  return nullptr;
}

template std::shared_ptr<P3MFFTBackend<float, P3MFFTKokkosConfig>>
make_p3m_kokkos_fft_backend<float>(boost::mpi::communicator const &,
                                   Utils::Vector3i const &,
                                   Utils::Vector3i const &,
                                   Utils::Vector3i const &,
                                   Utils::Vector3i const &);
template std::shared_ptr<P3MFFTBackend<double, P3MFFTKokkosConfig>>
make_p3m_kokkos_fft_backend<double>(boost::mpi::communicator const &,
                                    Utils::Vector3i const &,
                                    Utils::Vector3i const &,
                                    Utils::Vector3i const &,
                                    Utils::Vector3i const &);

#endif // ESPRESSO_P3M
