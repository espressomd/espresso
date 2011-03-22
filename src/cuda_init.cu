/*
  Copyright (C) 2010 The ESPResSo project
  
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
#include <cuda.h>

extern "C" {

#include "utils.h"

}

#include "cuda_init.h"

void gpu_init()
{
  int deviceCount, dev, found;
  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
    fprintf(stderr, "%d: cannot start CUDA.\n", this_node);
    errexit();
  }
  if (deviceCount == 0) {
    fprintf(stderr, "%d: There is no CUDA device.\n", this_node);
    errexit();
  }
  found = 0;
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major > 1 || (deviceProp.major == 1 && deviceProp.minor >= 1)) {
      fprintf(stderr, "%d:using CUDA device %s.\n", this_node, deviceProp.name);
      found = 1; break;
    }
  }
  if (!found) {
    fprintf(stderr, "%d: There is no device with compute capability >= 1.1.\n", this_node);
    errexit();
  }
}
