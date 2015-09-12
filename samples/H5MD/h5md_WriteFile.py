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
n_time=5
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
	#PARTICLE
	system.part[i].pos=numpy.array([0.1*i,0.2*i,0.3*i])
	system.part[i].v=numpy.array([0.2*i,0.3*i,0.4*i])
	system.part[i].f=numpy.array([0.3*i,0.4*i,0.5*i])
	system.part[i].type=i%3
	#system.part[i].bonds=[[outBond,(i+1)%n_part],[outBond,(i+2)%n_part]] #TODO
	system.part[i].mass=0.4*i
	system.part[i].omega_lab=numpy.array([0.5*i,0.6*i,0.7*i])
	system.part[i].rinertia=numpy.array([0.6*i,0.7*i,0.8*i])
	system.part[i].omega_body=numpy.array([0.7*i,0.8*i,0.9*i])
	system.part[i].torque_lab=numpy.array([0.8*i,0.9*i,1.0*i])
	system.part[i].quat=numpy.array([0.9*i,0.9*i,1.0*i,1.1*i])	
	####system.part[i].director=numpy.array([0.1,0.1,0.1,0.1])*[i*3,i*3,i*3,i*3]	#Not implemented yet
	system.part[i].q=1.0*i
	system.part[i].virtual=11*i          															#Only virtual or vs_relative in myconfig
	#system.part[i].vs_relative=numpy.array([1.2*i,1.3*i,1.4*i])	  #Only virtual or vs_relative in myconfig #ERROR in python code
	system.part[i].dip=numpy.array([1.3*i,1.4*i,1.5*i])
	system.part[i].dipm=14*i	
	system.part[i].ext_force=numpy.array([1.5*i,1.6*i,1.7*i])
	system.part[i].fix=numpy.array([i%2,i%2,i%2])
	system.part[i].ext_torque=numpy.array([1.6*i,1.7*i,1.8*i])
	system.part[i].gamma=1.7*i
	system.part[i].temp=1.8*i
# 	system.part[i].rotation=i%2	#Use only without analyze pressure
	####system.part[i].exclude=-1 #TODO
	####system.part[i].swimming=-1 #TODO
	system.box_l = numpy.array([100,200,300])
	#####system.part[i].image=numpy.array([0.1,0.1,0.1])*[i*2,i*2,i*2] #TOASK
	system.part[i].id=19*i	# Will be overwritten from Espresso


for i in range(n_time):
	system.time=0.1*(i+1)
	h5.write_to_h5.time(i,"particles/atoms/position/","time")
	h5.write_to_h5.time_step(i,"particles/atoms/position/","step")
	h5.write_to_h5.type(-1)
	h5.write_to_h5.pos(-1)
	h5.write_to_h5.v(-1)
	h5.write_to_h5.f(-1)
	#h5.write_to_h5.bonds(-1)
	h5.write_to_h5.mass(-1)
	h5.write_to_h5.omega_lab(-1)
	h5.write_to_h5.rinertia(-1)
	h5.write_to_h5.omega_body(-1)
	h5.write_to_h5.torque_lab(-1)
	h5.write_to_h5.quat(-1)
	####h5.write_to_h5.director(-1) #Not implemented yet
	h5.write_to_h5.q(-1)
	h5.write_to_h5.virtual(-1)      #Only virtual or vs_relative in myconfig
	#h5.write_to_h5.vs_relative(-1) #Only virtual or vs_relative in myconfig #ERROR in python code
	h5.write_to_h5.dip(-1)
	h5.write_to_h5.dipm(-1)
	h5.write_to_h5.ext_force(-1)
	h5.write_to_h5.fix(-1)
	h5.write_to_h5.ext_torque(-1)
	h5.write_to_h5.gamma(-1)
	h5.write_to_h5.temp(-1)
	h5.write_to_h5.rotation(-1)		#Use only without analyze pressure
	####h5.write_to_h5.exclude(-1) #TODO
	####h5.write_to_h5.swimming(-1) #TODO
	h5.write_to_h5.box_edges(-1)
	h5.write_to_h5.box_boundary(-1)
	h5.write_to_h5.box_dimension(-1)
	h5.write_to_h5.id(-1)
	#h5.write_to_h5.image(-1) #TOASK	
	h5.write_to_h5.userdefined(-1,[1,2,3],"User/user1/","value",'f8')
	
	#ANALYZE
	h5.write_to_h5.mindist(-1,[0],[1])
	h5.write_to_h5.distto(-1,None,[1,1,1])	
	h5.write_to_h5.analyze_linear_momentum(-1,True,True)
	h5.write_to_h5.nbhood(-1,[1,1,1],1.0,'3d')
	h5.write_to_h5.cylindrical_average(-1,[1,1,1],[1,1,1],1,1,1,1,[-1])
	h5.write_to_h5.pressure(-1,False)
	h5.write_to_h5.stress_tensor(-1,False)
	h5.write_to_h5.analyze_energy(-1,'all','default','default')
	h5.write_to_h5.calc_re(-1,1,1,1)
	h5.write_to_h5.calc_rg(-1,1,1,1)
	h5.write_to_h5.calc_rh(-1,1,1,1)
	h5.write_to_h5.structure_factor(-1,1,1)


h5.close()


# #H5MD: Write VMD parameters
# h5.h5_write_vmd_parameters()
# h5.h5_write_vmd_parameters_extra.chain(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.name(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.resid([1,2,3,4,5])
# h5.h5_write_vmd_parameters_extra.resname(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.segid(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.type(["A", "B", "C"])


#vmd
#test resize n_part+time
#bond time dependent
#Tab 4er
# image
#Test name ("A")
#Test read error meldung
# todo xxx 
#Manual
#option fur zeitunabhangig ---> -1

#check if h5md builds --> ifdef
#comment all
#save without time and test READ




####ERRORS
#torque_lab,quat,omega_lab,dipole: schreibt/liest beim h5 aufruf falsche werte
#vs_relative: vs_relative needs six args ---> vs_relative needs three args ; q = x[3] ---> q = x[2] ???


####TOASK
#mass time dependent
#energy etc only numbers
#box tensor
#bleibt analyzeargs.energy as usual