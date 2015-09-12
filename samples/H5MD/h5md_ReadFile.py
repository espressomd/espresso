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

from __future__ import print_function
import espressomd._system as es
import espressomd
from espressomd import integrate
from espressomd import h5md
from espressomd.interactions import HarmonicBond
import numpy
from __builtin__ import tuple

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
	system.part[i].pos=numpy.array([0,0,0])
	



h5.read_from_h5.time(n_time-1,"particles/atoms/position/","time")
h5.read_from_h5.type()
h5.read_from_h5.pos()
h5.read_from_h5.v()
h5.read_from_h5.f()
#h5.read_from_h5.bonds()
h5.read_from_h5.mass()
h5.read_from_h5.omega_lab()
h5.read_from_h5.rinertia()
h5.read_from_h5.omega_body()
h5.read_from_h5.torque_lab()
h5.read_from_h5.quat()
####h5.read_from_h5.director() #Not implemented yet
h5.read_from_h5.q()
h5.read_from_h5.virtual()      #Only virtual or vs_relative in myconfig
#h5.read_from_h5.vs_relative() #Only virtual or vs_relative in myconfig #ERROR in python code
h5.read_from_h5.dip()
h5.read_from_h5.dipm()
h5.read_from_h5.ext_force()
h5.read_from_h5.fix()
h5.read_from_h5.ext_torque()
h5.read_from_h5.gamma()
h5.read_from_h5.temp()
h5.read_from_h5.rotation()
####h5.read_from_h5.exclude() #TODO
####h5.read_from_h5.swimming() #TODO
h5.read_from_h5.box_edges()
h5.read_from_h5.id()
#h5.read_from_h5.image() #TOASK
result_user = h5.read_from_h5.userdefined(n_time-1,"User/user1/","value")








for i in range(n_part):
# 	print(system.time)
# 	print(system.part[i].type)
# 	print(system.part[i].pos)
# 	print(system.part[i].v)
# 	print(system.part[i].f)
# 	#print(system.part[i].bonds)
# 	print(system.part[i].mass)
# 	print(system.part[i].omega_lab)
# 	print(system.part[i].rinertia)
# 	print(system.part[i].omega_body)
# 	print(system.part[i].torque_lab)
# 	print(system.part[i].quat)
# 	####print(system.part[i].director)
# 	print(system.part[i].q)
# 	print(system.part[i].virtual)         															#Only virtual or vs_relative in myconfig
# 	#print(system.part[i].vs_relative)																  #Only virtual or vs_relative in myconfig #ERROR in python code
# 	print(system.part[i].dip)
# 	print(system.part[i].dipm)
#  	
# 	print(system.part[i].ext_force)
# 	print(system.part[i].fix)
# 	print(system.part[i].ext_torque)
# 	print(system.part[i].gamma)
# 	print(system.part[i].temp)
# 	print(system.part[i].rotation)
# 	####print(system.part[i].exclude) #TODO
# 	####print(system.part[i].swimming) #TODO
#  	
# 	print(system.box_l)
# 	#####print(system.part[i].image) #TOASK
# 	print(system.part[i].id)	
	print(result_user)













# #H5MD: Write VMD parameters
# h5.h5_write_vmd_parameters()
# h5.h5_write_vmd_parameters_extra.chain(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.name(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.resid([1,2,3,4,5])
# h5.h5_write_vmd_parameters_extra.resname(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.segid(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.type(["A", "B", "C"])


# h5.h5_read_vmd_parameters("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.chain("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.name("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.resid("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.resname("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.segid("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.type("A")#xxxxxxxxxxxxxxxx


#TODO def MASS  ,  remove es aus constructor  ,  remove reload file 
#TODO enegrie string to values
#TODO image
#attributes
#box_l nur [10,0,0]
