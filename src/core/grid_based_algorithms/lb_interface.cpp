/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

ActiveLB lattice_switch = ActiveLB::NONE;

namespace LB {

ActiveLB get_lattice_switch() { return lattice_switch; }

int get_steps_per_md_step(double md_timestep) {
  return static_cast<int>(std::round(get_tau() / md_timestep));
}

void init() {}

void propagate() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    lb_walberla()->integrate();
#endif
  }
}

void sanity_checks(double time_step) {
  if (lattice_switch == ActiveLB::NONE)
    return;

  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    lb_sanity_checks(*lb_walberla(), *lb_walberla_params(), time_step);
#endif
  }
}

void lebc_sanity_checks(unsigned int shear_direction,
                        unsigned int shear_plane_normal) {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    lb_walberla()->check_lebc(shear_direction, shear_plane_normal);
#endif
  }
}

double get_agrid() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    return lb_walberla_params()->get_agrid();
#endif
  }
  throw NoLBActive();
}

void check_tau_time_step_consistency(double tau, double time_step) {
  // use float epsilon since tau may be a float
  auto const eps = static_cast<double>(std::numeric_limits<float>::epsilon());
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

double get_tau() {
#ifdef WALBERLA
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    return lb_walberla_params()->get_tau();
  }
#endif
  throw NoLBActive();
}

double get_kT() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    return lb_walberla()->get_kT();
#endif
  }
  throw NoLBActive();
}

double get_lattice_speed() { return get_agrid() / get_tau(); }

#ifdef WALBERLA
namespace Walberla {

static Utils::Vector3d get_momentum() { return lb_walberla()->get_momentum(); }

static boost::optional<Utils::Vector3d>
get_velocity_at_pos(Utils::Vector3d pos) {
  return lb_walberla()->get_velocity_at_pos(pos);
}

REGISTER_CALLBACK_ONE_RANK(get_velocity_at_pos)

static boost::optional<double>
get_interpolated_density_at_pos(Utils::Vector3d pos) {
  return lb_walberla()->get_interpolated_density_at_pos(pos);
}

REGISTER_CALLBACK_ONE_RANK(get_interpolated_density_at_pos)

static Utils::VectorXd<9> get_pressure_tensor() {
  return lb_walberla()->get_pressure_tensor();
}

REGISTER_CALLBACK_REDUCTION(get_pressure_tensor, std::plus<>())

} // namespace Walberla
#endif // WALBERLA

Utils::VectorXd<9> const get_pressure_tensor() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    return ::Communication::mpiCallbacks().call(
        ::Communication::Result::reduction, std::plus<>(),
        Walberla::get_pressure_tensor);
#endif
  }
  throw NoLBActive();
}

Utils::Vector3d calc_fluid_momentum() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    return Walberla::get_momentum();
#endif
  }
  throw NoLBActive();
}

Utils::Vector3d const get_interpolated_velocity(Utils::Vector3d const &pos) {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_velocity_at_pos, folded_pos / get_agrid());
#endif
  }
  throw NoLBActive();
}

double get_interpolated_density(Utils::Vector3d const &pos) {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
#ifdef WALBERLA
    auto const folded_pos = folded_position(pos, box_geo);
    return mpi_call(::Communication::Result::one_rank,
                    Walberla::get_interpolated_density_at_pos,
                    folded_pos / get_agrid());
#endif
  }
  throw NoLBActive();
}

} // namespace LB
