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

#include "config.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "DipolarBarnesHut.hpp"
#include "../forces.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"
#include "energy.hpp"

#ifdef DIPOLAR_BARNES_HUT

void activate_dipolar_barnes_hut(float epssq, float itolsq)
{
    delete dipolarBarnesHut;
    dipolarBarnesHut = nullptr;
    // also necessary on 1 CPU or GPU, does more than just broadcasting
    mpi_bcast_coulomb_params();
    dipolarBarnesHut = new DipolarBarnesHut(espressoSystemInterface, epssq, itolsq);
    forceActors.push_back(dipolarBarnesHut);
    energyActors.push_back(dipolarBarnesHut);

    coulomb.Dmethod = DIPOLAR_BH_GPU;
}

void deactivate_dipolar_barnes_hut()
{
    if (dipolarBarnesHut)
    {
        forceActors.remove(dipolarBarnesHut);
        energyActors.remove(dipolarBarnesHut);
    }
        delete dipolarBarnesHut;
        dipolarBarnesHut = nullptr;
}

DipolarBarnesHut *dipolarBarnesHut = nullptr;

#endif

