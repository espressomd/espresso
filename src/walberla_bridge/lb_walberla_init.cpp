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
#include "LBWalberlaImpl.hpp"

#include "core/mpi/Environment.h"

#include <utils/Vector.hpp>

void walberla_mpi_init() {
  int argc = 0;
  char **argv = nullptr;
  static walberla::mpi::Environment m_env =
      walberla::mpi::Environment(argc, argv);
}

LBWalberlaBase *new_lb_walberla(
    const std::shared_ptr<walberla::WalberlaBlockForest> &blockforest,
    double viscosity, double density, double kT, unsigned int seed) {

  return new walberla::LBWalberlaImpl(blockforest, viscosity, density, 1u, kT,
                                      seed);
}
