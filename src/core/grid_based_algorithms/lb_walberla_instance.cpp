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
#include <lb_walberla_init.hpp>

#include <utils/Vector.hpp>

#include <memory>

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

void init_lb_walberla(double viscosity, double density, double agrid,
                      double tau, const Utils::Vector3d &box_dimensions,
                      const Utils::Vector3i &node_grid, double kT,
                      unsigned int seed) {
  // Exceptions need to be converted to runtime errors so they can be
  // handled from Python in a parallel simulation
  try {
    lb_walberla_instance = new_lb_walberla(viscosity, density, agrid, tau,
                                           box_dimensions, node_grid, kT, seed);
    lb_walberla_params_instance = new LBWalberlaParams{agrid, tau};
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
    lb_walberla_instance = nullptr;
    lb_walberla_params_instance = nullptr;
  }
}
REGISTER_CALLBACK(init_lb_walberla)

void destruct_lb_walberla() {
  delete lb_walberla_instance;
  delete lb_walberla_params_instance;
  lb_walberla_instance = nullptr;
  lb_walberla_params_instance = nullptr;
}
REGISTER_CALLBACK(destruct_lb_walberla)

void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau, double kT, unsigned int seed) {
  Communication::mpiCallbacks().call_all(init_lb_walberla, viscosity, density,
                                         agrid, tau, box_geo.length(),
                                         node_grid, kT, seed);
  if (lb_walberla_instance) {
    lb_lbfluid_set_lattice_switch(ActiveLB::WALBERLA);
    lb_lbfluid_sanity_checks();
  }
}

void mpi_destruct_lb_walberla() {
  lb_lbfluid_set_lattice_switch(ActiveLB::NONE);
  Communication::mpiCallbacks().call_all(destruct_lb_walberla);
}
#endif
