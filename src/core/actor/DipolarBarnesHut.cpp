/*
  Copyright (C) 2016,2017,2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "DipolarBarnesHut.hpp"
#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "forces.hpp"
#include "grid.hpp"

#ifdef DIPOLAR_BARNES_HUT

std::unique_ptr<DipolarBarnesHut> dipolarBarnesHut;

void activate_dipolar_barnes_hut(float epssq, float itolsq) {
  // also necessary on 1 CPU or GPU, does more than just broadcasting
  dipole.method = DIPOLAR_BH_GPU;
  mpi_bcast_coulomb_params();

  dipolarBarnesHut = std::make_unique<DipolarBarnesHut>(espressoSystemInterface,
                                                        epssq, itolsq);
  forceActors.push_back(dipolarBarnesHut.get());
  energyActors.push_back(dipolarBarnesHut.get());
}

void deactivate_dipolar_barnes_hut() {
  if (dipolarBarnesHut) {
    forceActors.remove(dipolarBarnesHut.get());
    energyActors.remove(dipolarBarnesHut.get());
    dipolarBarnesHut.reset();
  }
}

#endif
