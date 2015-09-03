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
system.time_step = 0.01
system.skin      = 0.4
n_part=10
n_time=10
int_steps=10

h5=h5md.h5md("h5md_WriteReadFile.h5",system)


###############################################
bondId=0
bondClass=HarmonicBond
params={"r_0":1.1, "k":5.2}
		
system.bondedInter[bondId]=bondClass(**params)
outBond=system.bondedInter[bondId]
tnIn=bondClass(**params).typeNumber()
tnOut=outBond.typeNumber()
outParams=outBond.params

#system.part[0].bonds=[[outBond,1],[outBond,2]]
################################################

for i in range(n_part):
	system.part[i].pos=numpy.array([0.1*i,0.1*i,0.1*i])
	system.part[i].v=numpy.array([0.2*i,0.2*i,0.2*i])
	system.part[i].f=numpy.array([0.3*i,0.3*i,0.3*i])
	system.part[i].type=i%3
	#system.part[i].bonds=[[outBond,(i+1)%n_part],[outBond,(i+2)%n_part]]
	system.part[i].mass=0.4*i
	system.part[i].omega_lab=numpy.array([0.5*i,0.5*i,0.5*i])
	system.part[i].rinertia=numpy.array([0.6*i,0.6*i,0.6*i])
	system.part[i].omega_body=numpy.array([0.7*i,0.7*i,0.7*i])
	system.part[i].torque_lab=numpy.array([0.8*i,0.8*i,0.8*i])
	system.part[i].quat=numpy.array([0.9*i,0.9*i,0.9*i,0.9*i])	
	####system.part[i].director=numpy.array([0.1,0.1,0.1,0.1])*[i*3,i*3,i*3,i*3]	#Not implemented yet
	system.part[i].q=1.0*i
	system.part[i].virtual=11*i          															#Only virtual or vs_relative in myconfig
	#system.part[i].vs_relative=numpy.array([1.2*i,1.2*i,1.2*i])	  #Only virtual or vs_relative in myconfig #ERROR in python code
	system.part[i].dip=numpy.array([1.3*i,1.3*i,1.3*i])
	system.part[i].dipm=14*i	
	
	system.part[i].ext_force=numpy.array([1.5*i,1.5*i,1.5*i])
	system.part[i].fix=numpy.array([i%2,i%2,i%2])
	system.part[i].ext_torque=numpy.array([1.6*i,1.6*i,1.6*i])
	system.part[i].gamma=1.7*i
	system.part[i].temp=1.8*i
	system.part[i].rotation=i%2
	####system.part[i].exclude=-1 #TODO
	####system.part[i].swimming=-1 #TODO
	
	#############################
	
	system.box_l = numpy.array([100,200,300])
	
	#####system.part[i].image=numpy.array([0.1,0.1,0.1])*[i*2,i*2,i*2] #TOASK
	#system.part[i].id=i*7

	







for i in range(n_time):
	system.time=0.1*(i+1)
	h5.h5_write_particles.type()
	h5.h5_write_particles.pos(i)
	h5.h5_write_particles.v(i)
	h5.h5_write_particles.f(i)
	#h5.h5_write_particles.bonds(i)
	h5.h5_write_particles.mass()
	h5.h5_write_particles.omega_lab(i)
	h5.h5_write_particles.rinertia(i)
	h5.h5_write_particles.omega_body(i)
	h5.h5_write_particles.torque_lab(i)
	h5.h5_write_particles.quat(i)
	####h5.h5_write_particles.director(i) #Not implemented yet
	h5.h5_write_particles.q(i)
	h5.h5_write_particles.virtual(i)     #Only virtual or vs_relative in myconfig
	#h5.h5_write_particles.vs_relative(i) #Only virtual or vs_relative in myconfig #ERROR in python code
	h5.h5_write_particles.dip(i)
	h5.h5_write_particles.dipm(i)
	h5.h5_write_particles.ext_force(i)
	h5.h5_write_particles.fix(i)
	h5.h5_write_particles.ext_torque(i)
	h5.h5_write_particles.gamma(i)
	h5.h5_write_particles.temp(i)
	h5.h5_write_particles.rotation(i)
	####h5.h5_write_particles.exclude(i) #TODO
	####h5.h5_write_particles.swimming(i) #TODO
	
	#####################################
	
	h5.h5_write_particles.box(i)
	#h5.h5_write_particles.image(i) #TOASK
	#h5.h5_write_observable(i,[2,3,4],"energies","Obs1")
	#h5.h5_write_particles.id()










# #H5MD: Write VMD parameters
# h5.h5_write_vmd_parameters()
# h5.h5_write_vmd_parameters_extra.chain(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.name(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.resid([1,2,3,4,5])
# h5.h5_write_vmd_parameters_extra.resname(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.segid(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.type(["A", "B", "C"])


#bond time dependent
#chunked
#Tab 4er
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
#Manual
#reihenfolge particle_data
#change case=time to -1,-1 ---> -1,time
#option fur zeitunabhangig ---> -1

####ERRORS
#torque_lab,quat: schreibt beim h5 aufruf falsche werte
#vs_relative: vs_relative needs six args ---> vs_relative needs three args ; q = x[3] ---> q = x[2] ???


####TOASK
#mass time dependent
