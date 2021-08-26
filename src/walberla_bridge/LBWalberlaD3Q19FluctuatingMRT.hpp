/*
 * Copyright (C) 2020 The ESPResSo project
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
#include "LBWalberlaImpl.hpp"
#ifdef __AVX2__
#include "generated_kernels/CollideSweepThermalizedAVX.h"
#define CollisionModelName walberla::pystencils::CollideSweepThermalizedAVX
#include "generated_kernels/FluctuatingMRTLatticeModelAvx.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModelAvx
#else
#include "generated_kernels/CollideSweepThermalized.h"
#define CollisionModelName walberla::pystencils::CollideSweepThermalized
#include "generated_kernels/FluctuatingMRTLatticeModel.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModel
#endif

namespace walberla {
class LBWalberlaD3Q19FluctuatingMRT
    : public LBWalberlaImpl<LatticeModelName, CollisionModelName> {
  using LatticeModel = LatticeModelName;

public:
  LBWalberlaD3Q19FluctuatingMRT(
      double viscosity, double density, const Utils::Vector3i &grid_dimensions,
      const Utils::Vector3i &node_grid, int n_ghost_layers, double kT,
      unsigned int seed,
      boost::optional<LeesEdwardsCallbacks> &&lees_edwards_callbacks)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers,
                       kT, seed, std::move(lees_edwards_callbacks)) {
    m_lattice_model = std::make_shared<LatticeModel>(
        m_last_applied_force_field_id, real_c(kT), -1., -1., -1., -1., seed, 0);
    setup_with_valid_lattice_model(density, seed, 0u);
  };
};

} // namespace walberla
#undef LatticeModelName
#undef CollisionModelName
