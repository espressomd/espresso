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
h5=h5md.h5md("h5md_WriteReadFile.h5",system)

n_part=10
n_time=10
int_steps=10


system.part[0].pos=numpy.array([0.2,0.2,0.2])
system.part[1].pos=numpy.array([0.3,0.3,0.3])

bondId=0
bondClass=HarmonicBond
params={"r_0":1.1, "k":5.2}
		
system.bondedInter[bondId]=bondClass(**params)
outBond=system.bondedInter[bondId]
tnIn=bondClass(**params).typeNumber()
tnOut=outBond.typeNumber()
outParams=outBond.params

print("XXX",outBond,tnIn,tnOut,outParams)

ttt=tuple([tuple([outBond,1]),tuple([outBond,2])])
#ttt=tuple([2,3])
for test in ttt:
	print(test)
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")


system.part[0].bonds=[[outBond,1],[outBond,2]]
#system.part[0].bonds=tuple([tuple([outBond,1]),tuple([outBond,2])])





for i in range(n_part):
	integrate.integrate(int_steps)
	system.part[i].pos=numpy.array([0.2,0.2,0.2])*[i*1,i*1,i*1]
	#####system.part[i].image=numpy.array([0.1,0.1,0.1])*[i*2,i*2,i*2]
	system.part[i].v=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	system.part[i].f=numpy.array([0.1,0.1,0.1])*[i*4,i*4,i*4]
	system.part[i].bonds=[0,1]
	#system.part[i].type=i*5
	#system.part[i].mass=i*6
	#system.part[i].omega_lab=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	#####system.part[i].rinertia=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	#system.part[i].omega_body=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	#system.part[i].torque_lab=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	#system.part[i].quat=numpy.array([0.1,0.1,0.1,0.1])*[i*3,i*3,i*3,i*3]
	#system.part[i].id=i*7
	#system.part[i].q=7
	#####system.part[i].virtual=2
	#system.part[i].vs_relative=[2,2.3]
	#system.part[i].dip=numpy.array([0.1,0.1,0.1])*[i*3,i*3,i*3]
	#system.part[i].dipm=i*7.1
system.box_l = [10,10,10]



for i in range(n_time):
	system.time=0.1*(i+1)
	h5.h5_write_particles.pos(i)
	#h5.h5_write_particles.image(i)
	h5.h5_write_particles.v(i)
	h5.h5_write_particles.f(i)
	#h5.h5_write_particles.omega_lab(i)
	######h5.h5_write_particles.rinertia(i)
	#h5.h5_write_particles.omega_body(i)
	#h5.h5_write_particles.torque_lab(i)
	#h5.h5_write_particles.quat(i)
	#h5.h5_write_particles.dip(i)
	#h5.h5_write_particles.dipm(i)
	#h5.h5_write_particles.box(i)
	#h5.h5_write_observable(i,[2,3,4],"energies","Obs1")
#h5.h5_write_particles.type()
#h5.h5_write_particles.mass()
#h5.h5_write_particles.id()
#h5.h5_write_particles.q()
#####h5.h5_write_particles.virtual()
#h5.h5_write_particles.vs_relative() #Only virtual or vs_relative in myconfig






# #H5MD: Write VMD parameters
# h5.h5_write_vmd_parameters()
# h5.h5_write_vmd_parameters_extra.chain(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.name(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.resid([1,2,3,4,5])
# h5.h5_write_vmd_parameters_extra.resname(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.segid(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.type(["A", "B", "C"])




#TODO def MASS  ,  remove es aus constructor  ,  remove reload file 
#TODO enegrie string to values
#TODO image
#attributes
#box_l nur [10,0,0]
#Test name ("A")
#Test read error meldung
# todo xxx 
#self_h5md_class
#bond
#box (x,1,1)(x,2,2)(x,3,3)
#get box dimension, tensor like
#self_h5md_class
#deaktivate MASS