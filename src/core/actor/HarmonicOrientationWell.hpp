/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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
#ifndef _ACTOR_HARMONICORIENTATIONWELL_HPP
#define _ACTOR_HARMONICORIENTATIONWELL_HPP

#include "config.hpp"

#ifdef CUDA
#ifdef ROTATION

#include "Actor.hpp"
#include "SystemInterface.hpp"
#include <iostream>

void HarmonicOrientationWell_kernel_wrapper(float x, float y, float z, float k, int n, float *quatu, float *torque);

class HarmonicOrientationWell : public Actor {
public:
  HarmonicOrientationWell(float x1, float x2, float x3, float _k, SystemInterface &s);

  virtual void computeTorques(SystemInterface &s) {
    HarmonicOrientationWell_kernel_wrapper(x,y,z,k,s.npart_gpu(), s.quatuGpuBegin(), s.torqueGpuBegin());
  };

  virtual ~HarmonicOrientationWell() {}
protected:
  float x,y,z;
  float k;
};

#endif
#endif
#endif
