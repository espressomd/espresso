/*
 * Copyright (C) 2018-2020 The ESPResSo project
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
#include "LBBoundary.hpp"

#include "communication.hpp"
#include "config.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <functional>
#include <stdexcept>

namespace LBBoundaries {
Utils::Vector3d LBBoundary::get_force() const {
#ifdef LB_BOUNDARIES
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const agrid = lb_lbfluid_get_agrid();
    Utils::Vector3d force{0, 0, 0};
    for (auto index_and_pos : lb_walberla()->node_indices_positions(true)) {
      // Convert to MD units
      auto const index = index_and_pos.first;
      auto const pos = index_and_pos.second * agrid;
      if (shape().is_inside(pos)) {
        auto node_force_density = lb_walberla()->get_node_boundary_force(index);
        if (node_force_density) {
          force += (*node_force_density);
        }
      }
    } // loop over lb cells
    return boost::mpi::all_reduce(comm_cart, force,
                                  std::plus<Utils::Vector3d>());
#endif
  }
#endif
  throw std::runtime_error("LB Boundary code called with inactive LB");
}
} // namespace LBBoundaries
