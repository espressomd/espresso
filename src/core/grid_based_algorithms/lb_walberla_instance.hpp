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

#include "config.hpp"

#ifdef LB_WALBERLA
#include "LbWalberlaBase.hpp"

/** @brief Initialize Walberla's MPI manager */
void walberla_mpi_init();

/** @brief Access the per-MPI-node LbWalberla isntance */
LbWalberlaBase *lb_walberla();

/** @brief Create the LbWalberla instance and sets the lattice switch to
 *  WALBERLA
 *
 *  @param viscosity Fluid viscosity
 *  @param density Fluid density
 *  @param agrid  Size of one lb cell
 *  @param tau    LB time step
 */
void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau, double kT, unsigned int seed);

/** @brief Destruct the LbWalberla instance and set lattice switch to NONE */
void mpi_destruct_lb_walberla();

#endif

#endif
