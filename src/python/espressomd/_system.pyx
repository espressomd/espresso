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


      

#
#      
#  def __init__(self, id):
#    self.id=id
#   
