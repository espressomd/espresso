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
#ifndef _HARMONICPOTENTIAL_HPP
#define _HARMONICPOTENTIAL_HPP

#include "config.hpp"

#ifdef CUDA

#include "Potential.hpp"
#include "EspressoSystemInterface.hpp"
#include <iostream>

void HarmonicPotential_kernel_wrapper(float x, float y, float z, float k,
		     int n, float *pos, float *f);

class HarmonicPotential : public Potential {
public:
  HarmonicPotential(float x1, float x2, float x3, float _k, SystemInterface &s) {
    x = x1;
    y = x2;
    z = x3;
    k = _k;

    if(!s.requestFGpu())
      std::cerr << "HarmonicPotential needs access to forces on GPU!" << std::endl;

    if(!s.requestRGpu())
      std::cerr << "HarmonicPotential needs access to positions on GPU!" << std::endl;

  }; 
  virtual void computeForces(SystemInterface &s) {
    HarmonicPotential_kernel_wrapper(x,y,z,k,s.npart_gpu(),
					 s.rGpuBegin(), s.fGpuBegin());
  };
protected:
  float x,y,z;
  float k;
};

inline void addHarmonicPotential(float x1, float x2, float x3, float _k) {
	potentials.push_back(new HarmonicPotential(x1, x2, x3, _k, espressoSystemInterface));
}

#endif
#endif
