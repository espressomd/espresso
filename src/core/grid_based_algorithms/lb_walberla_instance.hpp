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
#ifndef GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#define GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP

#include "config/config.hpp"

#ifdef WALBERLA
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include <memory>

struct LBWalberlaParams {
  LBWalberlaParams(double agrid, double tau) : m_agrid(agrid), m_tau(tau) {}
  double get_agrid() const { return m_agrid; };
  double get_tau() const { return m_tau; };

private:
  double m_agrid;
  double m_tau;
};

/** @brief Access the per-MPI-node waLBerla LB instance */
std::shared_ptr<LBWalberlaBase> lb_walberla();

/** @brief Access the waLBerla parameters */
std::shared_ptr<LBWalberlaParams> lb_walberla_params();

/** @brief Create the LBWalberla instance.
 *
 *  @param lb_lattice  LB lattice
 *  @param lb_params   LB parameters
 *  @param viscosity   Fluid viscosity
 *  @param density     Fluid density
 *  @param kT          Temperature
 *  @param seed        LB random seed
 *  @param single_precision   Use single precision
 */
std::shared_ptr<LBWalberlaBase>
init_lb_walberla(std::shared_ptr<LatticeWalberla> const &lb_lattice,
                 LBWalberlaParams const &lb_params, double viscosity,
                 double density, double kT, int seed, bool single_precision);

void lb_sanity_checks(LBWalberlaBase const &lb_fluid,
                      LBWalberlaParams const &lb_params, double md_time_step);

/** @brief Register a waLBerla LB instance and update lattice switch. */
bool activate_lb_walberla(std::shared_ptr<LBWalberlaBase> lb_fluid,
                          std::shared_ptr<LBWalberlaParams> lb_params);

/** @brief De-register a waLBerla LB instance and update lattice switch. */
void deactivate_lb_walberla();

#endif // WALBERLA

#endif
