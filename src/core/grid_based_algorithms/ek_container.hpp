/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_EK_CONTAINER_HPP
#define ESPRESSO_EK_CONTAINER_HPP

#include "config.hpp"

#ifdef LB_WALBERLA
#include "EKContainer.hpp"
#include "electrokinetics/EKinWalberlaBase.hpp"
#endif // LB_WALBERLA

namespace EK {
#ifdef LB_WALBERLA
extern EKContainer<EKinWalberlaBase> ek_container;
#endif // LB_WALBERLA

double get_tau();
int get_steps_per_md_step(double md_timestep);
void propagate();
} // namespace EK
#endif // ESPRESSO_EK_CONTAINER_HPP
