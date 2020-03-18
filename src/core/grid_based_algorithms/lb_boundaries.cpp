/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group,
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
/** \file
 *
 * Boundary conditions for lattice Boltzmann fluid dynamics.
 * Header file for \ref lb_boundaries.hpp.
 *
 */

#include "grid_based_algorithms/lb_boundaries.hpp"

#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "lbboundaries/LBBoundary.hpp"

#include <utils/index.hpp>
using Utils::get_linear_index;
#include <utils/constants.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

namespace LBBoundaries {

std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

void add(const std::shared_ptr<LBBoundary> &b) {
  lbboundaries.emplace_back(b);

  on_lbboundary_change();
}

void remove(const std::shared_ptr<LBBoundary> &b) {
  auto &lbb = lbboundaries;

  lbboundaries.erase(std::remove(lbb.begin(), lbb.end(), b), lbb.end());

  on_lbboundary_change();
}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
#if defined(LB_BOUNDARIES)

    lb_walberla()->clear_boundaries();

    auto const agrid = lb_lbfluid_get_agrid();

    for (auto index_and_pos : lb_walberla()->global_node_indices_positions()) {
      // Convert to MD units
      auto const index = index_and_pos.first;
      auto const pos = index_and_pos.second * agrid;

      for (auto it = lbboundaries.begin(); it != lbboundaries.end(); ++it) {
        double dist;
        Utils::Vector3d tmp;
        (**it).calc_dist(pos, dist, tmp);

        //        printf("%d %d %d, %g %g %g: %g %g %g:
        //        %g\n",index[0],index[1],index[2],pos[0],pos[1],pos[2],tmp[0],tmp[1],tmp[2],dist);
        if (dist <= 0) {

          // Set boundaries on the ghost layers
          auto const grid = lb_walberla()->get_grid_dimensions();
          for (int dx : {-1, 0, 1})
            for (int dy : {-1, 0, 1})
              for (int dz : {-1, 0, 1}) {
                Utils::Vector3i shifted_index =
                    index +
                    Utils::Vector3i{dx * grid[0], dy * grid[1], dz * grid[2]};
                lb_walberla()->set_node_velocity_at_boundary(
                    shifted_index,
                    (**it).velocity() / lb_lbfluid_get_lattice_speed());
              }
          break;
        } // if dist <=0
      }   // loop over boundaries
    }     // Loop over cells
#endif
#endif
  } // lattice switch is WALBERLA
}

#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */

} // namespace LBBoundaries
