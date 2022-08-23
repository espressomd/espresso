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

#include "../../include/walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp"

#include "../../include/walberla_bridge/LatticeWalberla.hpp"
#include "../../include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp"
#include "LBWalberlaImpl.hpp"

#include <core/mpi/Environment.h>

#include "utils/Vector.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

void walberla_mpi_init() {
  int argc = 0;
  char **argv = nullptr;
  static walberla::mpi::Environment m_env =
      walberla::mpi::Environment(argc, argv);
}

std::shared_ptr<LBWalberlaBase>
new_lb_walberla(std::shared_ptr<LatticeWalberla> const &lattice,
                double viscosity, double density, bool single_precision) {
  if (single_precision)
    return std::make_shared<walberla::LBWalberlaImpl<float>>(lattice, viscosity,
                                                             density);
  return std::make_shared<walberla::LBWalberlaImpl<double>>(lattice, viscosity,
                                                            density);
}

Utils::Vector3i calc_grid_dimensions(Utils::Vector3d const &box_size,
                                     double agrid) {
  Utils::Vector3i const grid_dimensions{
      static_cast<int>(std::round(box_size[0] / agrid)),
      static_cast<int>(std::round(box_size[1] / agrid)),
      static_cast<int>(std::round(box_size[2] / agrid))};
  for (int i : {0, 1, 2}) {
    if (std::abs(grid_dimensions[i] * agrid - box_size[i]) / box_size[i] >
        std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error(
          "Box length not commensurate with agrid in direction " +
          std::to_string(i) + " length " + std::to_string(box_size[i]) +
          " agrid " + std::to_string(agrid));
    }
  }
  return grid_dimensions;
}
