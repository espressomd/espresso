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

#include "LbWalberlaBase.hpp"
#include "LbWalberlaD3Q19FluctuatingMRT.hpp"
#include "LbWalberlaD3Q19MRT.hpp"

#include "core/mpi/Environment.h"

void walberla_mpi_init() {
  int argc = 0;
  char **argv = nullptr;
  static walberla::mpi::Environment m_env =
      walberla::mpi::Environment(argc, argv);
}

LbWalberlaBase *new_lb_walberla(double viscosity, double density, double agrid,
                                double tau,
                                const Utils::Vector3d &box_dimensions,
                                const Utils::Vector3i &node_grid, double kT,
                                unsigned int seed) {

  LbWalberlaBase *lb_walberla_instance;
  if (kT == 0.) { // un-thermalized LB
    lb_walberla_instance =
        new walberla::LbWalberlaD3Q19MRT(walberla::LbWalberlaD3Q19MRT{
            viscosity, density, agrid, tau, box_dimensions, node_grid, 1});
  } else { // thermalized LB
    lb_walberla_instance = new walberla::LbWalberlaD3Q19FluctuatingMRT(
        walberla::LbWalberlaD3Q19FluctuatingMRT{viscosity, density, agrid, tau,
                                                box_dimensions, node_grid, 1,
                                                kT, seed});
  }
  return lb_walberla_instance;
}
