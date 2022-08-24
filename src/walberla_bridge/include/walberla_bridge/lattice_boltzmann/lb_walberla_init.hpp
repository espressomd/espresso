/*
 * Copyright (C) 2019-2022 The ESPResSo project
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
#ifndef LB_WALBERLA_INIT_HPP
#define LB_WALBERLA_INIT_HPP

#include "LBWalberlaBase.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"

#include <utils/Vector.hpp>

#include <memory>

/** @brief Initialize Walberla's MPI manager */
void walberla_mpi_init();

std::shared_ptr<LBWalberlaBase>
new_lb_walberla(std::shared_ptr<LatticeWalberla> const &lattice,
                double viscosity, double density, bool single_precision);

Utils::Vector3i calc_grid_dimensions(Utils::Vector3d const &box_size,
                                     double agrid);

#endif
