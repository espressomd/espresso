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

#ifdef DIPOLAR_DIRECT_SUM

#include "DipolarDirectSum.hpp"

#include "EspressoSystemInterface.hpp"
#include "electrostatics_magnetostatics/common.hpp"
#include "energy.hpp"
#include "forces.hpp"

#include <memory>

static std::unique_ptr<DipolarDirectSum> dipolarDirectSum;

void activate_dipolar_direct_sum_gpu() {
  // also necessary on 1 CPU or GPU, does more than just broadcasting
  dipole.method = DIPOLAR_DS_GPU;
  mpi_bcast_coulomb_params();

  dipolarDirectSum =
      std::make_unique<DipolarDirectSum>(EspressoSystemInterface::Instance());
  forceActors.push_back(dipolarDirectSum.get());
  energyActors.push_back(dipolarDirectSum.get());
}

void deactivate_dipolar_direct_sum_gpu() {
  if (dipolarDirectSum) {
    forceActors.remove(dipolarDirectSum.get());
    energyActors.remove(dipolarDirectSum.get());
    dipolarDirectSum.reset();
  }
}

#endif
