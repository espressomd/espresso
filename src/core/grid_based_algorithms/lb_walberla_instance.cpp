/*
 * Copyright (C) 2019-2023 The ESPResSo project
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
#include "config/config.hpp"

#ifdef WALBERLA
#include "lb_walberla_instance.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>

#include <functional>
#include <memory>
#include <stdexcept>

static std::weak_ptr<LBWalberlaBase> lb_walberla_instance;
static std::shared_ptr<LBWalberlaParams> lb_walberla_params_instance;

std::shared_ptr<LBWalberlaBase> lb_walberla() {
  auto lb_walberla_instance_handle = ::lb_walberla_instance.lock();
  if (not lb_walberla_instance_handle) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberla instance.");
  }
  return lb_walberla_instance_handle;
}

std::shared_ptr<LBWalberlaParams> lb_walberla_params() {
  if (not ::lb_walberla_params_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberlaParams instance.");
  }
  return ::lb_walberla_params_instance;
}

void lb_sanity_checks(LBWalberlaBase const &lb_fluid,
                      LBWalberlaParams const &lb_params, double md_time_step) {
  auto const agrid = lb_params.get_agrid();
  auto const tau = lb_params.get_tau();
  // waLBerla and ESPResSo must agree on domain decomposition
  auto [lb_my_left, lb_my_right] = lb_fluid.get_lattice().get_local_domain();
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
    LB::check_tau_time_step_consistency(tau, md_time_step);
  }
}

void activate_lb_walberla(std::shared_ptr<LBWalberlaBase> lb_fluid,
                          std::shared_ptr<LBWalberlaParams> lb_params) {
  if (::lattice_switch != ActiveLB::NONE) {
    throw std::runtime_error("Cannot add a second LB instance");
  }
  lb_sanity_checks(*lb_fluid, *lb_params, get_time_step());
  auto const &lebc = ::box_geo.lees_edwards_bc();
  lb_fluid->check_lebc(lebc.shear_direction, lebc.shear_plane_normal);
  ::lb_walberla_instance = std::weak_ptr<LBWalberlaBase>{lb_fluid};
  ::lb_walberla_params_instance = lb_params;
  ::lattice_switch = ActiveLB::WALBERLA_LB;
}

void deactivate_lb_walberla() {
  ::lb_walberla_instance.reset();
  ::lb_walberla_params_instance.reset();
  ::lattice_switch = ActiveLB::NONE;
}

#endif
