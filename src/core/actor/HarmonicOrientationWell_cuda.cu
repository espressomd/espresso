/*
  Copyright (C) 2013,2014,2015,2016 The ESPResSo project
  
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

#include "cuda_utils.hpp"
#include <stdio.h>

__global__ void HarmonicOrientationWell_kernel(float x, float y, float z, float k,
				     int n, float *quatu, float *torque) {

  int id = blockIdx.x * blockDim.x + threadIdx.x;

  if(id >= n)
    return;

  float normdir = 1.0f/sqrt(x*x + y*y + z*z);
  float normori = 1.0f/sqrt(   quatu[3*id + 0]*quatu[3*id + 0]
                             + quatu[3*id + 1]*quatu[3*id + 1]
                             + quatu[3*id + 2]*quatu[3*id + 2] );

  float xn = x*normdir;
  float yn = y*normdir;
  float zn = z*normdir;

  float qn = quatu[3*id + 0]*normori;
  float rn = quatu[3*id + 1]*normori;
  float sn = quatu[3*id + 2]*normori;

  float sgndir = signbit(xn*qn + yn*rn + zn*sn);

  xn *= sgndir;
  yn *= sgndir;
  zn *= sgndir;

  torque[3*id + 0] += k*( -sn*yn + rn*zn );
  torque[3*id + 1] += k*( -qn*zn + sn*xn );
  torque[3*id + 2] += k*( -rn*xn + qn*yn );
}


void HarmonicOrientationWell_kernel_wrapper(float x, float y, float z, float k, int n, float *quatu, float *torque) {
  dim3 grid(1,1,1);
  dim3 block(1,1,1);

  if(n == 0)
    return;

  if(n <= 512) {
    grid.x = 1;
    block.x = n;
  } else {
    grid.x = n/512 + 1;
    block.x = 512;
  }

  KERNELCALL(HarmonicOrientationWell_kernel,grid,block,(x, y, z, k, n, quatu, torque))
}
