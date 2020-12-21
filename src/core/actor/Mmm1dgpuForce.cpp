/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifdef MMM1D_GPU

#include "actor/Mmm1dgpuForce.hpp"
#include "cells.hpp"
#include "energy.hpp"
#include "forces.hpp"
#include "grid.hpp"

#include <stdexcept>

void Mmm1dgpuForce::sanity_checks() {
  if (box_geo.periodic(0) || box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("MMM1D requires periodicity (0, 0, 1)");
  }
  if (cell_structure.decomposition_type() != CELL_STRUCTURE_NSQUARE) {
    throw std::runtime_error("MMM1D requires the N-square cellsystem");
  }
}

void Mmm1dgpuForce::activate() {
  forceActors.push_back(this);
  energyActors.push_back(this);
}

void Mmm1dgpuForce::deactivate() {
  forceActors.remove(this);
  energyActors.remove(this);
}

#endif
