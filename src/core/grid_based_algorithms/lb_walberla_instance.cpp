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

#include <LBWalberlaBase.hpp>
#include <LatticeWalberla.hpp>
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <cassert>
#include <functional>
#include <memory>
#include <stdexcept>

namespace {
std::shared_ptr<LBWalberlaBase> lb_walberla_instance = nullptr;
std::shared_ptr<LBWalberlaParams> lb_walberla_params_instance = nullptr;
} // namespace

std::shared_ptr<LBWalberlaBase> lb_walberla() {
  if (!lb_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberla instance.");
  }
  return lb_walberla_instance;
}

std::shared_ptr<LBWalberlaParams> lb_walberla_params() {
  if (!lb_walberla_params_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberlaParams instance.");
  }
  return lb_walberla_params_instance;
}

void mpi_activate_lb_walberla_local(
    std::shared_ptr<LBWalberlaBase> lb_fluid,
    std::shared_ptr<LBWalberlaParams> lb_params) {
  assert(::lattice_switch == ActiveLB::NONE);
  ::lb_walberla_instance = lb_fluid;
  ::lb_walberla_params_instance = lb_params;
  ::lattice_switch = ActiveLB::WALBERLA;
}

void mpi_deactivate_lb_walberla_local() {
  ::lb_walberla_instance = nullptr;
  ::lb_walberla_params_instance = nullptr;
  ::lattice_switch = ActiveLB::NONE;
}

static void lb_sanity_checks(LBWalberlaBase const &lb_fluid, double agrid) {
  // sanity check: waLBerla and ESPResSo must agree on domain decomposition
  auto [lb_my_left, lb_my_right] = lb_fluid.get_local_domain();
  lb_my_left *= agrid;
  lb_my_right *= agrid;
  auto const my_left = local_geo.my_left();
  auto const my_right = local_geo.my_right();
  auto const tol = agrid / 1E6;
  if ((lb_my_left - my_left).norm2() > tol or
      (lb_my_right - my_right).norm2() > tol) {
    runtimeErrorMsg() << this_node << ": left ESPResSo: [" << my_left << "], "
                      << "left waLBerla: [" << lb_my_left << "]\n";
    runtimeErrorMsg() << this_node << ": right ESPResSo: [" << my_right << "], "
                      << "right waLBerla: [" << lb_my_right << "]\n";
    throw std::runtime_error(
        "waLBerla and ESPResSo disagree about domain decomposition.");
  }
}

std::shared_ptr<LBWalberlaBase>
mpi_init_lb_walberla_local(std::shared_ptr<LatticeWalberla> lattice,
                           double viscosity, double density, double agrid,
                           double tau, double kT, unsigned int seed) {
  bool flag_failure = false;
  std::shared_ptr<LBWalberlaBase> lb_ptr;
  try {
    auto const md_time_step = get_time_step();
    if (md_time_step > 0.)
      check_tau_time_step_consistency(tau, md_time_step);
    lb_ptr = new_lb_walberla(lattice, viscosity, density, kT, seed);
    lb_sanity_checks(*lb_ptr, agrid);
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
