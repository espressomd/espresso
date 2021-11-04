/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include "ibm_volcons.hpp"

#include "communication.hpp"
#include "immersed_boundaries.hpp"

#include <cassert>
#include <stdexcept>

void mpi_set_n_ibm_volcons_bonds_local(int softID) {
  immersed_boundaries.register_softID(softID);
}

REGISTER_CALLBACK(mpi_set_n_ibm_volcons_bonds_local)

void mpi_set_n_ibm_volcons_bonds(int softID) {
  mpi_call_all(mpi_set_n_ibm_volcons_bonds_local, softID);
}

/** Set parameters of volume conservation */
IBMVolCons::IBMVolCons(const int softID, const double kappaV) {
  this->softID = softID;
  this->kappaV = kappaV;
  // NOTE: We cannot compute the reference volume here because not all
  // interactions are setup and thus we do not know which triangles belong to
  // this softID. Calculate it later in the init function of
  // \ref ImmersedBoundaries::init_volume_conservation()
  volRef = 0.;
  if (this_node == 0)
    mpi_set_n_ibm_volcons_bonds(softID);
}
