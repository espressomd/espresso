/*
 * Copyright (C) 2019-2020 The ESPResSo project
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

#include "lb_walberla_init.hpp"

#include "LBWalberlaBase.hpp"
#include "LBWalberlaD3Q19FluctuatingMRT.hpp"
#include "LBWalberlaD3Q19MRT.hpp"

#include "core/mpi/Environment.h"

#include <utils/Vector.hpp>

void walberla_mpi_init() {
  int argc = 0;
  char **argv = nullptr;
  static walberla::mpi::Environment m_env =
      walberla::mpi::Environment(argc, argv);
}

/*
#include "generated_kernels/MRTLatticeModel.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "lbm/field/Adaptors.h"
 */
LBWalberlaBase *new_lb_walberla(double viscosity, double density,
                                const Utils::Vector3i &grid_dimensions,
                                const Utils::Vector3i &node_grid, double kT,
                                unsigned int seed) {

  LBWalberlaBase *lb_walberla_instance;
  if (kT == 0.) { // un-thermalized LB
     auto ptr = new walberla::LBWalberlaD3Q19MRT(
        viscosity, density, grid_dimensions, node_grid, 1, kT, seed);
     lb_walberla_instance = ptr;
     /*
      * auto [bc, pdf_field] = ptr->foo({1,1,1});
     //< walberla::lbm::MRTLatticeModel>
     walberla::field::GhostLayerFieldAdaptor<walberla::lbm::DensityAdaptionFunction<walberla::lbm::MRTLatticeModel>, 0>().get( pdf_field->latticeModel(), *pdf_field, bc->cell.x(), bc->cell.y(), bc->cell.z());
      */
  } else { // thermalized LB
    lb_walberla_instance = new walberla::LBWalberlaD3Q19FluctuatingMRT(
        viscosity, density, grid_dimensions, node_grid, 1, kT, seed);
  }
  return lb_walberla_instance;
}
