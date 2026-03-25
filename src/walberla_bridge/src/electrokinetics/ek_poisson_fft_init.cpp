/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "PoissonSolverFFT.hpp"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

#include <waLBerlaDefinitions.h>

#if __has_include(<heffte.h>)
#define HAS_HEFFTE
#endif

#include <memory>
#include <stdexcept>

namespace walberla {

std::shared_ptr<walberla::PoissonSolver>
new_ek_poisson_fft(std::shared_ptr<LatticeWalberla> const &lattice,
                   double permittivity, bool single_precision) {
#if not defined(WALBERLA_BUILD_WITH_FFT) and not defined(HAS_HEFFTE)
  throw std::runtime_error("software was compiled without FFT support");
#else
  if (single_precision) {
    return std::make_shared<
        walberla::PoissonSolverFFT<float, lbmpy::Arch::CPU>>(lattice,
                                                             permittivity);
  }
  return std::make_shared<walberla::PoissonSolverFFT<double, lbmpy::Arch::CPU>>(
      lattice, permittivity);
#endif
}

} // namespace walberla
