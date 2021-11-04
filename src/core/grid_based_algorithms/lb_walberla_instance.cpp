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

#include <functional>
#include <memory>
#include <stdexcept>

namespace {
LBWalberlaBase *lb_walberla_instance = nullptr;
LBWalberlaParams *lb_walberla_params_instance = nullptr;
} // namespace

LBWalberlaBase *lb_walberla() {
  if (!lb_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberla instance.");
  }
  return lb_walberla_instance;
}

LBWalberlaParams *lb_walberla_params() {
  if (!lb_walberla_params_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LBWalberlaParams instance.");
  }
  return lb_walberla_params_instance;
}

void destruct_lb_walberla() {
  delete lb_walberla_instance;
  delete lb_walberla_params_instance;
  lb_walberla_instance = nullptr;
  lb_walberla_params_instance = nullptr;
}
REGISTER_CALLBACK(destruct_lb_walberla)

void mpi_init_lb_walberla_local(std::shared_ptr<LatticeWalberla> lattice,
                                double viscosity, double density, double agrid,
                                double tau, double kT, unsigned int seed) {
  bool flag_failure = false;
  // Exceptions need to be converted to runtime errors so they can be
  // handled from Python in a parallel simulation
  try {
    check_tau_time_step_consistency(tau, get_time_step());
    lb_walberla_instance =
        new_lb_walberla(lattice, viscosity, density, kT, seed);
    lb_walberla_params_instance = new LBWalberlaParams{agrid, tau};
    // TODO WALBERLA remove lattice_switch (modify lb_lbfluid_sanity_checks)
    ::lattice_switch = ActiveLB::WALBERLA;
    lb_lbfluid_sanity_checks(get_time_step());
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
    flag_failure = true;
  }
  boost::mpi::all_reduce(comm_cart, flag_failure, std::logical_or<>());
  if (flag_failure) {
    destruct_lb_walberla();
    // TODO WALBERLA remove lattice_switch
    ::lattice_switch = ActiveLB::NONE;
  }
}

void mpi_destruct_lb_walberla() {
  Communication::mpiCallbacks().call_all(destruct_lb_walberla);
}
#endif
