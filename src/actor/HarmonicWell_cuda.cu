/*
  Copyright (C) 2013,2014 The ESPResSo project
  
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
#include "HarmonicWell.hpp"

#include "cuda_utils.hpp"
#include <stdio.h>

__global__ void HarmonicWell_kernel(float x, float y, float z, float k,
				     int n, float *pos, float *f) {

  int id = blockIdx.x * blockDim.x + threadIdx.x;

  if(id >= n)
    return;

  f[3*id + 0] += k*(x - pos[3*id + 0]);
  f[3*id + 1] += k*(y - pos[3*id + 1]);
  f[3*id + 2] += k*(z - pos[3*id + 2]);
}


void HarmonicWell_kernel_wrapper(float x, float y, float z, float k, int n, float *pos, float *f) {
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

  KERNELCALL(HarmonicWell_kernel,grid,block,(x, y, z, k, n, pos, f))
}
