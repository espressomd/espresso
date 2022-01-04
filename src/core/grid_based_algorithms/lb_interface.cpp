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

#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
#include "grid_based_algorithms/lb_walberla_interface.hpp"

#include "BoxGeometry.hpp"
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace LB {
int get_steps_per_md_step(double md_timestep) {
  return static_cast<int>(std::round(lb_lbfluid_get_tau() / md_timestep));
}
} // namespace LB

ActiveLB lattice_switch = ActiveLB::NONE;

void lb_lbfluid_init() {}

void lb_lbfluid_propagate() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_walberla()->integrate();
#endif
  }
}

void lb_lbfluid_sanity_checks(double time_step) {
  if (lattice_switch == ActiveLB::NONE)
    return;

  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    lb_sanity_checks(*lb_walberla(), *lb_walberla_params(), time_step);
#endif
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla_params()->get_agrid();
#endif
  }
  throw NoLBActive();
}

void check_tau_time_step_consistency(double tau, double time_step) {
  auto const eps = std::numeric_limits<float>::epsilon();
  if ((tau - time_step) / (tau + time_step) < -eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be >= MD time_step (" +
                                std::to_string(time_step) + ")");
  auto const factor = tau / time_step;
  if (fabs(round(factor) - factor) / factor > eps)
    throw std::invalid_argument("LB tau (" + std::to_string(tau) +
                                ") must be an integer multiple of the "
                                "MD time_step (" +
                                std::to_string(time_step) + "). Factor is " +
                                std::to_string(factor));
}

double lb_lbfluid_get_tau() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    return lb_walberla_params()->get_tau();
  }
#endif
  throw NoLBActive();
}

double lb_lbfluid_get_kT() {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    return lb_walberla()->get_kT();
#endif
  }
  throw NoLBActive();
}

double lb_lbfluid_get_lattice_speed() {
  return lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
}

Utils::Vector6d lb_lbfluid_get_pressure_tensor_local() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    auto const &lb = lb_walberla();
    auto const gridsize = lb->lattice().get_grid_dimensions();
    Utils::Vector6d tensor{};
    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          const Utils::Vector3i node{{i, j, k}};
          auto const node_tensor = lb->get_node_pressure_tensor(node);
          if (node_tensor) {
            tensor += *node_tensor;
          }
        }
      }
    }
    return tensor;
  }
#endif
  return {};
}

REGISTER_CALLBACK_REDUCTION(lb_lbfluid_get_pressure_tensor_local,
                            std::plus<Utils::Vector6d>())

const Utils::Vector6d lb_lbfluid_get_pressure_tensor() {
#ifdef LB_WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA) {
    auto const visc = lb_walberla()->get_viscosity();
    auto const gridsize = lb_walberla()->lattice().get_grid_dimensions();
    auto const number_of_nodes = gridsize[0] * gridsize[1] * gridsize[2];
    auto tensor = ::Communication::mpiCallbacks().call(
        ::Communication::Result::reduction, std::plus<Utils::Vector6d>(),
        lb_lbfluid_get_pressure_tensor_local);
    tensor *= 1. / static_cast<double>(number_of_nodes);
    Walberla::walberla_off_diagonal_correction(tensor, visc);
    return tensor;
  }
#endif
  throw NoLBActive();
}

ActiveLB lb_lbfluid_get_lattice_switch() { return lattice_switch; }

Utils::Vector3d lb_lbfluid_calc_fluid_momentum() {
  Utils::Vector3d fluid_momentum{};
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    fluid_momentum = ::Communication::mpiCallbacks().call(
        ::Communication::Result::Reduction(), std::plus<>(),
        Walberla::get_momentum);
#endif
  } else
    throw NoLBActive();

  return fluid_momentum;
}

const Utils::Vector3d
lb_lbfluid_get_interpolated_velocity(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_velocity_at_pos,
                    folded_pos / lb_lbfluid_get_agrid());
#endif
  }
  throw NoLBActive();
}

double lb_lbfluid_get_interpolated_density(const Utils::Vector3d &pos) {
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_interpolated_density_at_pos,
                    folded_pos / lb_lbfluid_get_agrid());
#endif
  }
  throw NoLBActive();
}
