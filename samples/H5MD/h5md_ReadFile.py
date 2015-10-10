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

import espressomd._system as es
import espressomd
from espressomd import h5md
from espressomd.interactions import HarmonicBond
import numpy as np

# H5MD: Create/Read File dataset
#############################################################
system = espressomd.System()
h5_filename="h5md_WriteReadFile.h5"
h5=h5md.h5md("h5md_WriteReadFile.h5",system)

n_part=10
n_time=5
int_steps=10

#h5.h5_dataset_size(h5.h5_file['particles/position/step'])

for i in range(n_part):
    system.part[i].pos=np.array([0,0,0])

h5.read_from_h5.time(n_time-1,"particles/atoms/position/","time")
h5.read_from_h5.type(n_time-1)
h5.read_from_h5.pos(n_time-1)
h5.read_from_h5.v(n_time-1)
h5.read_from_h5.f(n_time-1)
#h5.read_from_h5.bonds(n_time-1)
h5.read_from_h5.mass(n_time-1)
h5.read_from_h5.omega_lab(n_time-1)
h5.read_from_h5.rinertia(n_time-1)
h5.read_from_h5.omega_body(n_time-1)
h5.read_from_h5.torque_lab(n_time-1)
h5.read_from_h5.quat(n_time-1)
####h5.read_from_h5.director(n_time-1) #Not implemented yet
h5.read_from_h5.q(n_time-1)
h5.read_from_h5.virtual(n_time-1)      #Only virtual or vs_relative in myconfig
#h5.read_from_h5.vs_relative(n_time-1) #Only virtual or vs_relative in myconfig #ERROR in python code
h5.read_from_h5.dip(n_time-1)
h5.read_from_h5.dipm(n_time-1)
h5.read_from_h5.ext_force(n_time-1)
h5.read_from_h5.fix(n_time-1)
h5.read_from_h5.ext_torque(n_time-1)
h5.read_from_h5.gamma(n_time-1)
h5.read_from_h5.temp(n_time-1)
h5.read_from_h5.rotation(n_time-1)
####h5.read_from_h5.exclude(n_time-1) #TODO
####h5.read_from_h5.swimming(n_time-1) #TODO
h5.read_from_h5.box_edges(n_time-1)
h5.read_from_h5.id(n_time-1)
#h5.read_from_h5.image(n_time-1) #TOASK
result_user = h5.read_from_h5.userdefined(n_time-1,"User/user1/","value")








for i in range(n_part):
#     print(system.time)
#     print(system.part[i].type)
#     print(system.part[i].pos)
    print(system.part[i].v)
#     print(system.part[i].f)
#     #print(system.part[i].bonds)
#     print(system.part[i].mass)
#     print(system.part[i].omega_lab)
#     print(system.part[i].rinertia)
#     print(system.part[i].omega_body)
#     print(system.part[i].torque_lab)
#     print(system.part[i].quat)
#     ####print(system.part[i].director)
#     print(system.part[i].q)
# 	  print(system.part[i].virtual)         															#Onl     virtual or vs_relative in myconfig
#     #print(system.part[i].vs_relative)																  #Only virtual or vs_relative in myconfig #ERROR in python code
#     print(system.part[i].dip)
#     print(system.part[i].dipm)
#      
#     print(system.part[i].ext_force)
#     print(system.part[i].fix)
#     print(system.part[i].ext_torque)
#     print(system.part[i].gamma)
#     print(system.part[i].temp)
#     print(system.part[i].rotation)
#     ####print(system.part[i].exclude) #TODO
#     ####print(system.part[i].swimming) #TODO
#      
#     print(system.box_l)
#     #####print(system.part[i].image) #TOASK
#     print(system.part[i].id)    
#     print(result_user)



