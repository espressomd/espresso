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
#include "config.hpp"

#ifdef LB_WALBERLA
#include "lb_walberla_instance.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"
#include "lees_edwards_protocol.hpp"

#include <LBWalberlaBase.hpp>
#include <LatticeWalberla.hpp>
#include <LeesEdwardsPack.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>

#include <cassert>
#include <functional>
#include <memory>
#include <stdexcept>

// TODO walberla: use proper script interface logic
std::weak_ptr<LeesEdwards::ActiveProtocol> lees_edwards_active_protocol{};

namespace {
static std::weak_ptr<LBWalberlaBase> lb_walberla_instance;
static std::shared_ptr<LBWalberlaParams> lb_walberla_params_instance;
} // namespace

std::shared_ptr<LBWalberlaBase> lb_walberla() {
  auto lb_walberla_instance_handle = lb_walberla_instance.lock();
  if (!lb_walberla_instance_handle) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberla instance.");
  }
  return lb_walberla_instance_handle;
}

std::shared_ptr<LBWalberlaParams> lb_walberla_params() {
  if (!lb_walberla_params_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberlaParams instance.");
  }
  return lb_walberla_params_instance;
}

void lb_sanity_checks(LBWalberlaBase const &lb_fluid,
                      LBWalberlaParams const &lb_params, double md_time_step) {
  auto const agrid = lb_params.get_agrid();
  auto const tau = lb_params.get_tau();
  // waLBerla and ESPResSo must agree on domain decomposition
  auto [lb_my_left, lb_my_right] = lb_fluid.lattice().get_local_domain();
  lb_my_left *= agrid;
  lb_my_right *= agrid;
  auto const my_left = local_geo.my_left();
  auto const my_right = local_geo.my_right();
  auto const tol = agrid / 1E6;
  if ((lb_my_left - my_left).norm2() > tol or
      (lb_my_right - my_right).norm2() > tol) {
    runtimeErrorMsg() << "\nMPI rank " << this_node << ": "
                      << "left ESPResSo: [" << my_left << "], "
                      << "left waLBerla: [" << lb_my_left << "]"
                      << "\nMPI rank " << this_node << ": "
                      << "right ESPResSo: [" << my_right << "], "
                      << "right waLBerla: [" << lb_my_right << "]";
    throw std::runtime_error(
        "waLBerla and ESPResSo disagree about domain decomposition.");
  }
  // LB time step and MD time step must agree
  if (md_time_step > 0.) {
    check_tau_time_step_consistency(tau, md_time_step);
  }
}

bool activate_lb_walberla(std::shared_ptr<LBWalberlaBase> lb_fluid,
                          std::shared_ptr<LBWalberlaParams> lb_params) {
  bool flag_failure = false;
  try {
    assert(::lattice_switch == ActiveLB::NONE);
    lb_sanity_checks(*lb_fluid, *lb_params, get_time_step());
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "during waLBerla activation: " << e.what();
    flag_failure = true;
  }
  if (boost::mpi::all_reduce(comm_cart, flag_failure, std::logical_or<>())) {
    return true;
  }
  ::lb_walberla_instance = std::weak_ptr<LBWalberlaBase>{lb_fluid};
  ::lb_walberla_params_instance = lb_params;
  ::lattice_switch = ActiveLB::WALBERLA;
  return false;
}

void deactivate_lb_walberla() {
  ::lb_walberla_instance.reset();
  ::lb_walberla_params_instance.reset();
  ::lattice_switch = ActiveLB::NONE;
}

std::shared_ptr<LBWalberlaBase>
init_lb_walberla(std::shared_ptr<LatticeWalberla> const &lb_lattice,
                 LBWalberlaParams const &lb_params, double viscosity,
                 double density, double kT, int seed, bool single_precision) {
  bool flag_failure = false;
  std::shared_ptr<LBWalberlaBase> lb_ptr;
  try {
    assert(seed >= 0);
    lb_ptr = new_lb_walberla(lb_lattice, viscosity, density, single_precision);
    if (kT != 0.)
      lb_ptr->set_collision_model(kT, seed);
    if (auto active_protocol = lees_edwards_active_protocol.lock()) {
      auto const &le_bc = box_geo.clees_edwards_bc();
      auto lees_edwards_object = LeesEdwardsPack(
          le_bc.shear_direction, le_bc.shear_plane_normal,
          [active_protocol]() {
            return get_pos_offset(get_sim_time(), *active_protocol);
          },
          [active_protocol]() {
            return get_shear_velocity(get_sim_time(), *active_protocol);
          });
      lb_ptr->set_collision_model(std::move(lees_edwards_object));
    }
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "during waLBerla initialization: " << e.what();
    flag_failure = true;
  }
  if (boost::mpi::all_reduce(comm_cart, flag_failure, std::logical_or<>())) {
    return nullptr;
  }
  return lb_ptr;
}
#endif
