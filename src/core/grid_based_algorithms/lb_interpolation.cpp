/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include <boost/mpi/collectives.hpp>

#include "communication.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include <utils/Vector.hpp>

namespace {

InterpolationOrder interpolation_order = InterpolationOrder::linear;
}

void mpi_set_interpolation_order(InterpolationOrder const &order) {
  interpolation_order = order;
}

REGISTER_CALLBACK(mpi_set_interpolation_order)

void lb_lbinterpolation_set_interpolation_order(
    InterpolationOrder const &order) {
  mpi_call_all(mpi_set_interpolation_order, order);
}

InterpolationOrder lb_lbinterpolation_get_interpolation_order() {
  return interpolation_order;
}

const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &pos) {
  Utils::Vector3d interpolated_u{};

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA

    auto res = lb_walberla()->get_velocity_at_pos(pos / lb_lbfluid_get_agrid());
    if (!res) {
      auto folded_pos = folded_position(pos, box_geo);
      res = lb_walberla()->get_velocity_at_pos(folded_pos /
                                               lb_lbfluid_get_agrid());
    }

    if (!res) {
      printf("%d: positoin: %g %g %g\n", this_node, pos[0], pos[1], pos[2]);
      throw std::runtime_error(
          "Interpolated velocity could not be obtained from Walberla");
    }
    extern double sim_time;
    return *res;
#endif
  }
  throw std::runtime_error("No LB active.");
}

void lb_lbinterpolation_add_force_density(
    const Utils::Vector3d &pos, const Utils::Vector3d &force_density) {
  switch (interpolation_order) {
  case (InterpolationOrder::quadratic):
    throw std::runtime_error(
        "The non-linear interpolation scheme is not implemented.");
  case (InterpolationOrder::linear):
    if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
      lb_walberla()->add_force_at_pos(pos / lb_lbfluid_get_agrid(),
                                      force_density);
#endif
    } else
      throw std::runtime_error("No LB active.");
  }
}
