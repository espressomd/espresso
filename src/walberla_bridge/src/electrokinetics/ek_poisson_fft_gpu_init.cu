/*
 * Copyright (C) 2025 The ESPResSo project
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

#if defined(__NVCC__)
#define RESTRICT __restrict__
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 554 // no implicit or explicit cast
#else
#pragma push
#pragma diag_suppress 554 // no implicit or explicit cast
#endif
#endif

#include "FFT_CUDA.cuh"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/PoissonSolver/PoissonSolverCuFFT.hpp>
#include <walberla_bridge/electrokinetics/ek_poisson_fft_gpu_init.hpp>

#include <memory>

namespace walberla {

std::shared_ptr<walberla::PoissonSolver>
new_ek_poisson_fft_cuda(std::shared_ptr<LatticeWalberla> const &lattice,
                        double permittivity, bool single_precision) {
  if (single_precision) {
    return std::make_shared<walberla::PoissonSolverCuFFT<float>>(lattice,
                                                                 permittivity);
  }
  return std::make_shared<walberla::PoissonSolverCuFFT<double>>(lattice,
                                                                permittivity);
}

} // namespace walberla
