/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  
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

#include "utils.hpp"
#include "parser.hpp"
#include "cuda_init.hpp"

#ifdef CUDA

/** prints a list of the available GPUs to the result of the Tcl interpreter.
    Only devices with compute capability of 1.1 or higher are listed, since
    atomic operations are required for CUDA-LB. */
static int list_gpus(Tcl_Interp *interp)
{
  int deviceCount = cuda_get_n_gpus();
  if (deviceCount < 0) {
    Tcl_AppendResult(interp, "cannot initialize CUDA: ", cuda_error, (char *)NULL);
    return TCL_ERROR;
  }

  int found = 0; 
  for (int dev = 0; dev < deviceCount; ++dev) {
    // look only for devices with compute capability > 1.1 (for atomic operations)
    if (cuda_check_gpu(dev) == ES_OK) {
      char id[4 + 64 + TCL_INTEGER_SPACE];
      char name[64];
      cuda_get_gpu_name(dev, name);
      sprintf(id, " {%d %.64s}", dev, name);
      Tcl_AppendResult(interp, id, NULL);
      found = 1;
    }
  }
  if (found == 0) {
    Tcl_AppendResult(interp, "no GPUs present", NULL);
  }

  return TCL_OK;
}
#endif /* defined(CUDA) */

/** returns 1 if and only if the GPU with the given id is usable for
    CUDA computations.  Only devices with compute capability of 1.1 or
    higher are ok, since atomic operations are required for
    CUDA-LB. */
int tclcommand_cuda(ClientData data, Tcl_Interp *interp,
		    int argc, char **argv)
{
#ifndef CUDA
    Tcl_AppendResult(interp, "Feature CUDA required!", (char *)NULL);
    return TCL_ERROR;
#else
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
    if (argc <= 1 || !ARG1_IS_I(dev)) {
      Tcl_AppendResult(interp, "expected: cuda setdevice <devnr>", (char *)NULL);
      return TCL_ERROR;
    }
    if (cuda_check_gpu(dev) == ES_ERROR) {
      Tcl_AppendResult(interp, "GPU not present or compute model not sufficient", (char *)NULL);
      return TCL_ERROR;
    }
    if (cuda_set_device(dev) == ES_OK) {
      return TCL_OK;
    }
    else {
      Tcl_AppendResult(interp, cuda_error, (char *)NULL);
      return TCL_ERROR;
    }
  }
  else if (ARG0_IS_S("getdevice")) {
    if (argc != 1) {
      Tcl_AppendResult(interp, "cuda getdevice takes no arguments", (char *)NULL);
      return TCL_ERROR;
    }
    int dev = cuda_get_device();
    if (dev >= 0) {
      char buffer[TCL_INTEGER_SPACE];
      sprintf(buffer, "%d", dev);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      return TCL_OK;
    }
    else {
      Tcl_AppendResult(interp, cuda_error, (char *)NULL);
      return TCL_ERROR;
    }
  }
  else {
    Tcl_AppendResult(interp, "unknown subcommand \"", argv[0], "\"", (char *)NULL);
    return TCL_ERROR;
  }
#endif /* defined(CUDA) */
}
