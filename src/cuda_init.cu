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

// CUDA code is always interpreted as C++, so we need the extern C interface
extern "C" {

#include "utils.h"
#include "parser.h"
#include "cuda_init.h"

}

static int list_gpus(Tcl_Interp *interp)
{
  int deviceCount, dev;

  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
    Tcl_AppendResult(interp, "cannot initialize CUDA", NULL);
    return TCL_ERROR;
  }

  // look for devices with compute capability > 1.1 (for atomic operations)
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major > 1 || (deviceProp.major == 1 && deviceProp.minor >= 1)) {
      char id[4 + 64 + TCL_INTEGER_SPACE];
      sprintf(id, " {%d %.64s}", dev, deviceProp.name);
      Tcl_AppendResult(interp, id, NULL);
    }
  }
  return TCL_OK;
}

int tclcommand_cuda(ClientData data, Tcl_Interp *interp,
		    int argc, char **argv)
{
  if (argc <= 1) {
    Tcl_AppendResult(interp, "too few arguments to the cuda command", (char *)NULL);
    return TCL_ERROR;
  }
  argc--; argv++;
  
  if (ARG0_IS_S("list")) {
    if (argc != 1) {
      Tcl_AppendResult(interp, "cuda list takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    return list_gpus(interp);
  }
  else if (ARG0_IS_S("setdevice")) {
    int dev;
    cudaError_t error;
    if (argc <= 1 || !ARG1_IS_I(dev)) {
      Tcl_AppendResult(interp, "expected: cuda setdevice <devnr>", (char *)NULL);
      return TCL_ERROR;
    }
    error = cudaSetDevice(dev);
    if (error == cudaSuccess) {
      return TCL_OK;
    }
    else {
      Tcl_AppendResult(interp, cudaGetErrorString(error), (char *)NULL);
      return TCL_ERROR;
    }
  }
  else if (ARG0_IS_S("getdevice")) {
    if (argc != 1) {
      Tcl_AppendResult(interp, "cuda getdevice takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    int dev;
    cudaError_t error;
    error = cudaGetDevice(&dev);
    if (error == cudaSuccess) {
      char buffer[TCL_INTEGER_SPACE];
      sprintf(buffer, "%d", dev);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      return TCL_OK;
    }
    else {
      Tcl_AppendResult(interp, cudaGetErrorString(error), (char *)NULL);
      return TCL_ERROR;
    }
  }
  else {
    Tcl_AppendResult(interp, "unknown subcommand \"", argv[0], "\"", (char *)NULL);
    return TCL_ERROR;
  }
}
