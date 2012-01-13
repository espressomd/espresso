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
#ifndef CUDA_INIT_H
#define CUDA_INIT_H

/** external functions called by the Tcl-command to set the CUDA device to use or retrieve information
    available devices. */
int check_gpu(int dev);
int list_gpus(Tcl_Interp *interp);
int setdevice(int dev);
int getdevice(int* dev);
int getdevicelist(int* devl, char* devname);
#endif
