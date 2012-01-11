/*
  Copyright (C) 2010,2011 The ESPResSo project
  
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
#include "utils.h"
#include "parser.h"
// CUDA code is always interpreted as C++, so we need the extern C interface
extern "C" {

#include "cuda_init.h"

}

/** \name minimally required compute capability. */
/*@{*/
static const int computeCapabilityMinMajor = 1;
static const int computeCapabilityMinMinor = 1;
/*@}*/

/** returns 1 if and only if the GPU with the given id is usable for
    CUDA computations.  Only devices with compute capability of 1.1 or
    higher are ok, since atomic operations are required for
    CUDA-LB. */
int check_gpu(int dev)
{
  cudaDeviceProp deviceProp;
  if (cudaGetDeviceProperties(&deviceProp, dev) != cudaSuccess) {
    return 0;
  }
  if (deviceProp.major < computeCapabilityMinMajor ||
      (deviceProp.major == computeCapabilityMinMajor &&
       deviceProp.minor < computeCapabilityMinMinor)) {
    return 0;
  }
  return 1;
}

/** prints a list of the available GPUs to the result of the Tcl interpreter.
    Only devices with compute capability of 1.1 or higher are listed, since
    atomic operations are required for CUDA-LB. */
int list_gpus(Tcl_Interp *interp)
{
  int deviceCount, dev;

  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
    Tcl_AppendResult(interp, "cannot initialize CUDA", NULL);
    return TCL_ERROR;
  }

  // look for devices with compute capability > 1.1 (for atomic operations)
  for (dev = 0; dev < deviceCount; ++dev) {
    if (check_gpu(dev)) {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      char id[4 + 64 + TCL_INTEGER_SPACE];
      sprintf(id, " {%d %.64s}", dev, deviceProp.name);
      Tcl_AppendResult(interp, id, NULL);
    }
  }
  return TCL_OK;
}

