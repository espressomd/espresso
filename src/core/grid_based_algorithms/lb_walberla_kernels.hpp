/*
 * Copyright (C) 2019-2022 The ESPResSo project
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
#ifndef LB_WALBERLA_KERNELS_HPP
#define LB_WALBERLA_KERNELS_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "grid_based_algorithms/lb_walberla_interface.hpp"

#include <LBWalberlaBase.hpp>

#include <utils/Vector.hpp>

namespace Walberla {

/**
 * @brief Compute the local contribution to the average pressure tensor.
 * Needs to be summed over all MPI ranks.
 */
inline Utils::Vector6d average_pressure_tensor_local(LBWalberlaBase const &lb) {
  auto const grid_size = lb.lattice().get_grid_dimensions();
  auto const number_of_nodes = Utils::product(grid_size);
  Utils::Vector6d tensor{};
  for (int i = 0; i < grid_size[0]; i++) {
    for (int j = 0; j < grid_size[1]; j++) {
      for (int k = 0; k < grid_size[2]; k++) {
        auto const index = Utils::Vector3i{{i, j, k}};
        auto const node_tensor = lb.get_node_pressure_tensor(index);
        if (node_tensor) {
          tensor += *node_tensor;
        }
      }
    }
  }
  tensor *= 1. / static_cast<double>(number_of_nodes);
  walberla_off_diagonal_correction(tensor, lb.get_viscosity());
  return tensor;
}

} // namespace Walberla
#endif
#endif
