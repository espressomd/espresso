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

import espressomd
from espressomd import h5md
import numpy as np

### Prepare system for h5-reading
system = espressomd.System()

### Prepare particles for h5-reading
n_part=10
for i in range(n_part):
    system.part.add(id=i, pos=np.array([0,0,0]))
    
#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------TIME DEPENDENT-----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
### Open h5-file
n_time=5
h5=h5md.h5md("File.h5",system)

### Read from h5-file
# user defined 
result_user1 = h5.read_from_h5.userdefined("User/user1/","value1",(3,))    # One dimensional (pay attention to the comma (3,) )
result_user2 = h5.read_from_h5.userdefined("User/user1/","value2",(3,4,5)) # Three dimensional

# Particle data 
h5.read_from_h5.pos(n_time-1)
#h5.read_from_h5.v(n_time-1)
#h5.read_from_h5.f(n_time-1)
#h5.read_from_h5.mass(n_time-1)
#h5.read_from_h5.omega_lab(n_time-1)
#h5.read_from_h5.rinertia(n_time-1)
#h5.read_from_h5.omega_body(n_time-1)
#h5.read_from_h5.torque_lab(n_time-1)
#h5.read_from_h5.quat(n_time-1)
#h5.read_from_h5.q(n_time-1)
#h5.read_from_h5.virtual(n_time-1)      
#h5.read_from_h5.vs_relative(n_time-1) 
#h5.read_from_h5.dip(n_time-1)
#h5.read_from_h5.dipm(n_time-1)
#h5.read_from_h5.ext_force(n_time-1)
#h5.read_from_h5.fix(n_time-1)
#h5.read_from_h5.ext_torque(n_time-1)
#h5.read_from_h5.gamma(n_time-1)
#h5.read_from_h5.temp(n_time-1)
#h5.read_from_h5.rotation(n_time-1)
#h5.read_from_h5.box_edges(n_time-1)
#h5.read_from_h5.id(n_time-1)
#h5.read_from_h5.time(n_time-1,"particles/atoms/position/","time")
#h5.read_from_h5.type()

### Output results
print("\nTIME DEPENDENT")
print(result_user1)
print(result_user2)  
for i in range(n_part):  
    print(system.part[i].pos)
    #print(system.part[i].v)
    #print(system.part[i].f)
    #print(system.part[i].mass)
    #print(system.part[i].omega_lab)
    #print(system.part[i].rinertia)
    #print(system.part[i].omega_body)
    #print(system.part[i].torque_lab)
    #print(system.part[i].quat)
    #print(system.part[i].q)
 	  #print(system.part[i].virtual)         															
    #print(system.part[i].vs_relative)																  
    #print(system.part[i].dip)
    #print(system.part[i].dipm)
    #print(system.part[i].ext_force)
    #print(system.part[i].fix)
    #print(system.part[i].ext_torque)
    #print(system.part[i].gamma)
    #print(system.part[i].temp)
    #print(system.part[i].rotation)
    #print(system.box_l)
    #print(system.part[i].id)  
    #print(system.time)  
    #print(system.part[i].type)
    


#----------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------TIME INDEPENDENT----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
### Open h5-file
h5_time_independent=h5md.h5md("File_time_independent.h5",system)

### Read from h5-file
h5_time_independent.read_from_h5.pos()
h5_time_independent.read_from_h5.v(-1,"particles/atoms/velocity/","value")    #With additional parameters, use "-1" as first parameter 

### Output results
print("\nTIME INDEPENDENT")
for i in range(n_part):
    print(system.part[i].pos)
    print(system.part[i].v)

