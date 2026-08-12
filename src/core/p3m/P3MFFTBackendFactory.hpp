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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "p3m/P3MFFTBackend.hpp"
#include "p3m/common.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <boost/mpi/communicator.hpp>

#include <memory>

/** @brief The only FFT configuration the kokkos-fft backend implements:
 *  row-major layout, r2c transform reducing the contiguous last axis.
 */
using P3MFFTKokkosConfig =
    P3MFFTConfig<Utils::MemoryOrder::ROW_MAJOR, Utils::MemoryOrder::ROW_MAJOR,
                 true, 2u>;

/**
 * @brief Build the single-rank kokkos-fft backend (@c P3MFFTKokkos) when it
 * can serve this run, i.e. when ESPResSo was built with kokkos-fft support and
 * the communicator spans a single MPI rank; return @c nullptr otherwise so the
 * caller falls back to heFFTe.
 *
 * Defined once, in a translation unit of the @c espresso_p3m target -- the
 * only scope where the @c ESPRESSO_KOKKOS_FFT compile definition exists. The
 * templated P3M solver calls this function instead of branching on the macro
 * itself, so every translation unit that instantiates the solver (including
 * unit tests, which are compiled without the macro) agrees on one definition
 * and on the backend-selection behaviour.
 */
template <typename FloatType>
std::shared_ptr<P3MFFTBackend<FloatType, P3MFFTKokkosConfig>>
make_p3m_kokkos_fft_backend(boost::mpi::communicator const &comm,
                            Utils::Vector3i const &global_mesh,
                            Utils::Vector3i const &rs_local_ld_index,
                            Utils::Vector3i const &rs_local_ur_index,
                            Utils::Vector3i const &node_grid);

#endif // ESPRESSO_P3M
