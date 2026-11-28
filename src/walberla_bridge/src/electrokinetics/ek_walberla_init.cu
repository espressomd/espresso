/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

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

#include "EKinWalberlaImpl.hpp"

#if defined(__NVCC__)
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic pop
#else
#pragma pop
#endif
#endif

#include <walberla_bridge/Architecture.hpp>

#include "reactions/EKReactionImplBulk.hpp"
#include "reactions/EKReactionImplIndexed.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBaseIndexed.hpp>

#include <utils/Vector.hpp>

#include <waLBerlaDefinitions.h>

#include <memory>
#include <stdexcept>

namespace walberla {

std::shared_ptr<EKinWalberlaBase>
new_ek_walberla_gpu(std::shared_ptr<LatticeWalberla> const &lattice,
                    double diffusion, double kT, double valency,
                    Utils::Vector3d ext_efield, double density, bool advection,
                    bool friction_coupling, bool single_precision,
                    bool thermalized, unsigned int seed) {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  throw std::runtime_error("waLBerla was compiled without CUDA support");
#else
  if (single_precision) {
    return std::make_shared<EKinWalberlaImpl<13, float, lbmpy::Arch::GPU>>(
        lattice, diffusion, kT, valency, ext_efield, density, advection,
        friction_coupling, thermalized, seed);
  }

  return std::make_shared<EKinWalberlaImpl<13, double, lbmpy::Arch::GPU>>(
      lattice, diffusion, kT, valency, ext_efield, density, advection,
      friction_coupling, thermalized, seed);
#endif
}

std::shared_ptr<EKReactionBase> new_ek_reaction_bulk_gpu(
    std::shared_ptr<LatticeWalberla> const &lattice,
    typename EKReactionBase::reactants_type const &reactants,
    double coefficient) {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  throw std::runtime_error("waLBerla was compiled without CUDA support");
#else
  return std::make_shared<EKReactionImplBulk<lbmpy::Arch::GPU>>(
      lattice, reactants, coefficient);
#endif
}

std::shared_ptr<EKReactionBaseIndexed> new_ek_reaction_indexed_gpu(
    std::shared_ptr<LatticeWalberla> const &lattice,
    typename EKReactionBase::reactants_type const &reactants,
    double coefficient) {
#if not defined(WALBERLA_BUILD_WITH_CUDA)
  throw std::runtime_error("waLBerla was compiled without CUDA support");
#else
  return std::make_shared<EKReactionImplIndexed<lbmpy::Arch::GPU>>(
      lattice, reactants, coefficient);
#endif
}

} // namespace walberla
