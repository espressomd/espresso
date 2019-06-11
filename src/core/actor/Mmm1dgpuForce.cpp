/*
Copyright (C) 2010-2018 The ESPResSo project

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
#include "actor/Mmm1dgpuForce.hpp"
#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "forces.hpp"
#include "grid.hpp"

#ifdef MMM1D_GPU

void Mmm1dgpuForce::check_periodicity() {
  if (PERIODIC(0) || PERIODIC(1) || !PERIODIC(2)) {
    throw std::runtime_error("MMM1D requires periodicity (0,0,1)");
  }
}

#endif
