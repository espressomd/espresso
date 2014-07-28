#
# Copyright (C) 2013,2014 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
include "myconfig.pxi"
cimport numpy as np
import numpy as np
import particle_data
cimport particle_data 
import interactions
cimport interactions
import global_variables
from integrate import integrate
import thermostat
cimport thermostat
from changeVolume import changeVolume
from invalidateSystem import invalidateSystem
import cellsystem
cimport cellsystem
import analyze
cimport analyze
import utils
cimport utils

import debye_hueckel
#import lb
cimport cuda_init
import cuda_init
#cimport myconfig

import code_info

#public enum:
#  ERROR=-1


### Now come all the includes

cdef int instance_counter 

### Here we make a minimalistic Tcl_Interp available
cdef extern from "tcl.h":
  cdef struct Tcl_Interp:
    pass
    
  Tcl_Interp* Tcl_CreateInterp()
  int Tcl_Eval(Tcl_Interp * interp, char* script)
  char* Tcl_GetStringResult(Tcl_Interp * interp)

## Let's import on_program_start
cdef extern from *:
  void mpi_init(int* argc, char*** argv)
  int this_node
  void mpi_stop()

cdef mpi_init_helper():
  cdef int i
  cdef char** c
  i=0
  c=NULL
  mpi_init(&i, &c)

cdef extern from "initialize.hpp":
  void on_program_start()
  void mpi_loop()

## Now comes the real deal
cdef class EspressoHandle:
  cdef Tcl_Interp* interp
  cdef public int this_node
  def __init__(self):
    global instance_counter
    if instance_counter >= 1:
      raise Exception("Espresso shall only be instanciated once!")
    else:
      instance_counter+=1
      mpi_init_helper()
      self.this_node=this_node
      if this_node==0:
        self.interp = Tcl_CreateInterp() 
        self.Tcl_Eval('global argv; set argv ""')
        self.Tcl_Eval('set tcl_interactive 0')
        on_program_start()
      else:
        on_program_start()
        mpi_loop()
        self.die()

  def __del__(self):
    self.die()
    raise Exception("Espresso can not be deleted")

  def Tcl_Eval(self, string):
    result = Tcl_Eval(self.interp, string)
    tclresult = Tcl_GetStringResult(self.interp)
    if result:
      raise Exception("Tcl reports an error", tclresult)
  
    return tclresult

  def die(self):
    mpi_stop()

_espressoHandle=EspressoHandle()
glob = global_variables.GlobalsHandle()
part = particle_data.particleList()
#lbfluid=lb.DeviceList()
IF CUDA == 1:
    cu=cuda_init.CudaInitHandle()

# def TclEval(string):
#   if instance_counter == 0:
#     raise Exception("Espresso not initialized")
#   if instance_counter == 1:
#     _espressoHandle.Tcl_Eval(string)
nonBondedInter = interactions.NonBondedInteractions()
bondedInter = interactions.BondedInteractions()

if this_node==0:
  glob = global_variables.GlobalsHandle()
else:
  #why exit all other
  exit()

      

#
#      
#  def __init__(self, id):
#    self.id=id
#   
