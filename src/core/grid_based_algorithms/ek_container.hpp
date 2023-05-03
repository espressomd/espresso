/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA
#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#endif // WALBERLA

#include <stdexcept>

struct NoEKActive : public std::exception {
  const char *what() const noexcept override { return "EK not activated"; }
};

namespace EK {

#ifdef WALBERLA
extern EKContainer<EKinWalberlaBase> ek_container;
#endif // WALBERLA

double get_tau();
int get_steps_per_md_step(double md_timestep);
void propagate();

} // namespace EK

#endif
