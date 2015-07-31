/*
  Copyright (C) 2014 The ESPResSo project

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
#include "HarmonicOrientationWell.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"

#ifdef CUDA
#ifdef ROTATION

HarmonicOrientationWell::
HarmonicOrientationWell(float x1, float x2, float x3, float _k, SystemInterface &s)
  : x(x1), y(x2), z(x3), k(_k){
  if(!s.requestQuatuGpu())
    std::cerr << "HarmonicOrientationWell needs access to quatu on GPU!" << std::endl;

  if(!s.requestTorqueGpu())
    std::cerr << "HarmonicOrientationWell needs access to torques on GPU!" << std::endl;
}

#endif
#endif
