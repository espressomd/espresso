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
#include "cuda_init_tcl.h"

}

/** returns 1 if and only if the GPU with the given id is usable for
    CUDA computations.  Only devices with compute capability of 1.1 or
    higher are ok, since atomic operations are required for
    CUDA-LB. */

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
    if (!check_gpu(dev)) {
      Tcl_AppendResult(interp, "GPU not present or compute model not sufficient", (char *)NULL);
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
  else if (ARG0_IS_S("listdevice")) {
    if (argc != 1) {
      Tcl_AppendResult(interp, "cuda listdevice takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    int gpu_n;
    cudaError_t error;
    error = cudaGetDeviceCount(&gpu_n);
    if (error == cudaSuccess) {
      char buffer[TCL_INTEGER_SPACE];
      sprintf(buffer, "%d", gpu_n);
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
