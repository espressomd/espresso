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
#ifndef LB_WALBERLA_INTERFACE_HPP
#define LB_WALBERLA_INTERFACE_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <LBWalberlaNodeState.hpp>

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <vector>

namespace Walberla {

void clear_boundaries();
void update_boundary_from_shape(std::vector<int> const &raster_flat,
                                std::vector<double> const &slip_velocity_flat);
void update_boundary_from_list(std::vector<int> const &nodes_flat,
                               std::vector<double> const &vel_flat);

void set_node_from_checkpoint(Utils::Vector3i ind, LBWalberlaNodeState cpt);
boost::optional<LBWalberlaNodeState> get_node_checkpoint(Utils::Vector3i ind);
void do_reallocate_ubb_field();
void do_ghost_communication();

Utils::Vector3d get_momentum();

boost::optional<Utils::Vector3d> get_velocity_at_pos(Utils::Vector3d pos);
boost::optional<double> get_interpolated_density_at_pos(Utils::Vector3d pos);

inline void walberla_off_diagonal_correction(Utils::Vector6d &tensor,
                                             double visc) {
  auto const revert_factor = visc / (visc + 1.0 / 6.0);
  tensor[1] *= revert_factor;
  tensor[3] *= revert_factor;
  tensor[4] *= revert_factor;
}

} // namespace Walberla
#endif
#endif
