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
 * Source file for \ref lb_boundaries.hpp.
 */

#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"

#include "event.hpp"
#include "lbboundaries/LBBoundary.hpp"
#include "walberla_blockforest.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <vector>

namespace LBBoundaries {

std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES)

void add(const std::shared_ptr<LBBoundary> &b) {
  auto &lbb = lbboundaries;
  assert(std::find(lbb.begin(), lbb.end(), b) == lbb.end());
  lbb.emplace_back(b);

  on_lbboundary_change();
}

void remove(const std::shared_ptr<LBBoundary> &b) {
  auto &lbb = lbboundaries;
  assert(std::find(lbb.begin(), lbb.end(), b) != lbb.end());
  lbb.erase(std::remove(lbb.begin(), lbb.end(), b), lbb.end());

  on_lbboundary_change();
}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_walberla()->clear_boundaries();

    auto const agrid = lb_lbfluid_get_agrid();

    for (auto index_and_pos : lb_walberla()->node_indices_positions(true)) {
      // Convert to MD units
      auto const index = index_and_pos.first;
      auto const pos = index_and_pos.second * agrid;

      for (auto const &lbboundary : lbboundaries) {
        if (lbboundary->shape().is_inside(pos)) {
          lb_walberla()->set_node_velocity_at_boundary(
              index, lbboundary->velocity() / lb_lbfluid_get_lattice_speed(),
              false);
        }
      }
    }
    lb_walberla()->reallocate_ubb_field();
#endif
  } // lattice switch is WALBERLA
}

#endif /* LB_BOUNDARIES */

} // namespace LBBoundaries
