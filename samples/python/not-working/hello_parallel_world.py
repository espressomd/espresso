#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from __future__ import print_function

import espresso as es
import numpy
import code_info

print(code_info.electrostatics_defined())

es._espressoHandle.Tcl_Eval("thermostat lb 1.")
dev="cpu"
N=100
es.glob.time_step=0.01
es.glob.skin=1.
es.glob.box_l=[10., 10., 10.]
#print es.cu.device_list
#es.cu.device=0
es.lbfluid[dev].agrid=1
es.lbfluid[dev].dens=1
es.lbfluid[dev].visc=1
es.lbfluid[dev].friction=1
es.lbfluid[dev].tau=0.1

es.lbfluid[dev].ext_force=[1., 2., 3.,]
print(es.lbfluid[dev].ext_force)
print(es.lbfluid[dev].dens)
print(es.lbfluid[dev].visc)
print(es.lbfluid[dev].agrid)
es.lbfluid[dev].print_vtk_velocity="test.vtk"
#es.lb[dev].checkpoint_style=1
#es.lb[dev].checkpoint="cp.dat"

for i in range(N):
  es.part[i].pos=numpy.random.random(3)*es.glob.box_l

es.inter[0,0].lennardJones = {"eps":1,"sigma":1,"shift":0.25}

es._espressoHandle.Tcl_Eval("integrate 100")
#for i in range(N):
#  print es.part[i].pos

es._espressoHandle.die()

