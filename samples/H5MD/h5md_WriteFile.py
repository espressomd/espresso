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
    system.part[i].pos=np.array([0.1*i,0.2*i,0.3*i])
    system.part[i].v=np.array([0.2*i,0.3*i,0.4*i])
    system.part[i].f=np.array([0.3*i,0.4*i,0.5*i])
    system.part[i].type=i%3
    #system.part[i].bonds=[[outBond,(i+1)%n_part],[outBond,(i+2)%n_part]] #TODO
    system.part[i].mass=0.4*i
    system.part[i].omega_lab=np.array([0.5*i,0.6*i,0.7*i])
    system.part[i].rinertia=np.array([0.6*i,0.7*i,0.8*i])
    system.part[i].omega_body=np.array([0.7*i,0.8*i,0.9*i])
    system.part[i].torque_lab=np.array([0.8*i,0.9*i,1.0*i])
    system.part[i].quat=np.array([0.9*i,0.9*i,1.0*i,1.1*i])
    ####system.part[i].director=np.array([0.1,0.1,0.1,0.1])*[i*3,i*3,i*3,i*3]    #Not implemented yet
    system.part[i].q=1.0*i
    system.part[i].virtual=11*i          															#Only virtual or vs_relative in myconfig
    #system.part[i].vs_relative=np.array([1.2*i,1.3*i,1.4*i])	  #Only virtual or vs_relative in myconfig #ERROR in python code
    system.part[i].dip=np.array([1.3*i,1.4*i,1.5*i])
    system.part[i].dipm=14*i    
    system.part[i].ext_force=np.array([1.5*i,1.6*i,1.7*i])
    system.part[i].fix=np.array([i%2,i%2,i%2])
    system.part[i].ext_torque=np.array([1.6*i,1.7*i,1.8*i])
    system.part[i].gamma=1.7*i
    system.part[i].temp=1.8*i
#     system.part[i].rotation=i%2	#Use only without analyze pressure
    ####system.part[i].exclude=-1 #TODO
    ####system.part[i].swimming=-1 #TODO
    system.box_l = np.array([100,200,300])
    #####system.part[i].image=np.array([0.1,0.1,0.1])*[i*2,i*2,i*2] #TOASK
    system.part[i].id=19*i    # Will be overwritten from Espresso
  
for i in range(n_time):
    system.time=0.1*(i+1)
    h5.write_to_h5.time(i,"particles/atoms/position/","time")
    h5.write_to_h5.time_step(i,"particles/atoms/position/","step")
    h5.write_to_h5.type(i)
    h5.write_to_h5.pos(i)
    h5.write_to_h5.v(i)
    h5.write_to_h5.f(i)
    #h5.write_to_h5.bonds(i)
    h5.write_to_h5.mass(i)
    h5.write_to_h5.omega_lab(i)
    h5.write_to_h5.rinertia(i)
    h5.write_to_h5.omega_body(i)
    h5.write_to_h5.torque_lab(i)
    h5.write_to_h5.quat(i)
    ####h5.write_to_h5.director(i) #Not implemented yet
    h5.write_to_h5.q(i)
    h5.write_to_h5.virtual(i)      #Only virtual or vs_relative in myconfig
    #h5.write_to_h5.vs_relative(i) #Only virtual or vs_relative in myconfig #ERROR in python code
    h5.write_to_h5.dip(i)
    h5.write_to_h5.dipm(i)
    h5.write_to_h5.ext_force(i)
    h5.write_to_h5.fix(i)
    h5.write_to_h5.ext_torque(i)
    h5.write_to_h5.gamma(i)
    h5.write_to_h5.temp(i)
    h5.write_to_h5.rotation(i)        #Use only without analyze pressure
    ####h5.write_to_h5.exclude(i) #TODO
    ####h5.write_to_h5.swimming(i) #TODO
    h5.write_to_h5.box_edges(i)
    h5.write_to_h5.box_boundary(i)
    h5.write_to_h5.box_dimension(i)
    h5.write_to_h5.id(i)
    #h5.write_to_h5.image(i) #TOASK    
    h5.write_to_h5.userdefined(i,[1,2,3],"User/user1/","value",'f8')
    
    #ANALYZE
    h5.write_to_h5.mindist(i,[0],[1])
    h5.write_to_h5.distto(i,None,[1,1,1])    
    h5.write_to_h5.analyze_linear_momentum(i,True,True)
    h5.write_to_h5.nbhood(i,[1,1,1],1.0,'3d')
    h5.write_to_h5.cylindrical_average(i,[1,1,1],[1,1,1],1,1,1,1,[-1])
    h5.write_to_h5.pressure(i,False)
    h5.write_to_h5.stress_tensor(i,False)
    h5.write_to_h5.analyze_energy(i,'all','default','default')
    h5.write_to_h5.calc_re(i,1,1,1)
    h5.write_to_h5.calc_rg(i,1,1,1)
    h5.write_to_h5.calc_rh(i,1,1,1)
    h5.write_to_h5.structure_factor(i,1,1)
h5.close()


#vmd
# image,bonds
# todo xxx 
#Manual
#check if h5md builds --> ifdef
#comment all


####ERRORS
#torque_lab,quat,omega_lab,dipole: schreibt/liest beim h5 aufruf falsche werte
#vs_relative: vs_relative needs six args ---> vs_relative needs three args ; q = x[3] ---> q = x[2] ???


####TOASK
#box tensor
#bleibt analyzeargs.energy as usual