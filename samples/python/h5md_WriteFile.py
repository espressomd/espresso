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
from espressomd.interactions import HarmonicBond
import numpy as np

# Prepare system for h5-writing
system = espressomd.System()
system.time_step = 0.01
system.skin      = 0.4
n_part=10
n_time=5
int_steps=10

### Prepare bonds for h5-writing
bondId=0
bondClass=HarmonicBond
params={"r_0":1.1, "k":5.2}
system.bonded_inter[bondId]=bondClass(**params)
outBond=system.bonded_inter[bondId]
outBond = system.bonded_inter[bondId]
tnIn = bondClass(**params).type_number()
tnOut = outBond.type_number()
outParams = outBond.params

### Prepare particles for h5-writing
for i in range(n_part):
    system.part.add(id=i, pos=np.array([0,0,0]))
    
#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------TIME DEPENDENT-----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
### Create h5-file
h5=h5md.h5md("File.h5",system)

### Write to h5-file
for j in range(n_time):
    ### Set arbitrary values for writing
    system.time=0.1*(j+1)
    for i in range(n_part):
        system.part[i].pos=np.array([0.1*(i+j),0.2*(i+j),0.3*(i+j)])
        #system.part[i].v=np.array([0.2*(i+j),0.3*(i+j),0.4*(i+j)])
        #system.part[i].f=np.array([0.3*(i+j),0.4*(i+j),0.5*(i+j)])
        #system.part[i].bonds=[[outBond,1],[outBond,2],[outBond,3]]
        #system.part[i].mass=0.1*(i+j)
        #system.part[i].omega_lab=np.array([0.5*(i+j),0.6*(i+j),0.7*(i+j)])
        #system.part[i].rinertia=np.array([0.6*(i+j),0.7*(i+j),0.8*(i+j)])
        #system.part[i].omega_body=np.array([0.7*(i+j),0.8*(i+j),0.9*(i+j)])
        #system.part[i].torque_lab=np.array([0.8*(i+j),0.9*(i+j),1.0*(i+j)])
        #system.part[i].quat=np.array([0.9*(i+j),0.9*(i+j),1.0*(i+j),1.1*(i+j)])
        #system.part[i].q=1.0*(i+j)
        #system.part[i].virtual=11*(i+j)                                        
        #system.part[i].vs_relative=np.array([1.2*(i+j),1.3*(i+j),1.4*(i+j)])    
        #system.part[i].dip=np.array([1.3*(i+j),1.4*(i+j),1.5*(i+j)])
        #system.part[i].dipm=14*(i+j)    
        #system.part[i].ext_force=np.array([1.5*(i+j),1.6*(i+j),1.7*(i+j)])
        #system.part[i].fix=np.array([i%2,i%2,i%2])
        #system.part[i].ext_torque=np.array([1.6*(i+j),1.7*(i+j),1.8*(i+j)])
        #system.part[i].gamma=1.7*(i+j)
        #system.part[i].temp=1.8*(i+j)
        #system.part[i].rotation=i%2  
        #system.box_l = np.array([100,200,300])
        #system.part[i].type=i%3
      
    ### Write to h5-file 
    # user defined 
    h5.write_to_h5.userdefined(123.4,"User/user1/","value1",'f8',(3,))  # One dimensional (pay attention to the comma (3,) )
    h5.write_to_h5.userdefined(567.8,"User/user1/","value2",'f8',(3,4,5),(1,1,5)) # Three dimensional and manual chunk size
    
    # Particle data 
    h5.write_to_h5.pos(j)
    #h5.write_to_h5.v(j)
    #h5.write_to_h5.f(j)
    #h5.write_to_h5.bonds(j)
    #h5.write_to_h5.mass(j)
    #h5.write_to_h5.omega_lab(j)
    #h5.write_to_h5.rinertia(j)
    #h5.write_to_h5.omega_body(j)
    #h5.write_to_h5.torque_lab(j)
    #h5.write_to_h5.quat(j)
    #h5.write_to_h5.q(j)
    #h5.write_to_h5.virtual(j)      
    #h5.write_to_h5.vs_relative(j) 
    #h5.write_to_h5.dip(j)
    #h5.write_to_h5.dipm(j)
    #h5.write_to_h5.ext_force(j)
    #h5.write_to_h5.fix(j)
    #h5.write_to_h5.ext_torque(j)
    #h5.write_to_h5.gamma(j)
    #h5.write_to_h5.temp(j)
    #h5.write_to_h5.rotation(j)        
    #h5.write_to_h5.box_edges(j)
    #h5.write_to_h5.box_boundary(j)
    #h5.write_to_h5.box_dimension(j)
    #h5.write_to_h5.id(j)
    #h5.write_to_h5.time(j,"particles/atoms/position/","time")
    #h5.write_to_h5.time_step(j,"particles/atoms/position/","step")
    #h5.write_to_h5.type()
    
    # Observables 
    h5.write_to_h5.analyze_energy(j,'all','default','default')
    #h5.write_to_h5.mindist(j,[0],[1])
    #h5.write_to_h5.distto(j,None,[1,1,1])    
    #h5.write_to_h5.analyze_linear_momentum(j,True,True)
    #h5.write_to_h5.nbhood(j,[1,1,1],1.0,'3d')
    #h5.write_to_h5.cylindrical_average(j,[1,1,1],[1,1,1],1,1,1,1,[-1])
    #h5.write_to_h5.pressure(j,False)
    #h5.write_to_h5.stress_tensor(j,False)   
    #h5.write_to_h5.calc_re(j,1,1,1)
    #h5.write_to_h5.calc_rg(j,1,1,1)
    #h5.write_to_h5.calc_rh(j,1,1,1)
    #h5.write_to_h5.structure_factor(j,1,1)
    
# VMD 
h5.write_to_h5.VMD('indexOfSpecies',[0,1,2])
#h5.write_to_h5.VMD('species')
#h5.write_to_h5.VMD('charge')
#h5.write_to_h5.VMD('mass')
#h5.write_to_h5.VMD('bond_from')
#h5.write_to_h5.VMD('bond_to')
#h5.write_to_h5.VMD('radius',[1.0,2.0,3.0])
#h5.write_to_h5.VMD('name',['H','He','Li'])
#h5.write_to_h5.VMD('type',['H','He','Li'])
#h5.write_to_h5.VMD('resname',['r1','r2'])
#h5.write_to_h5.VMD('resid',[0,0,0,0,0,1,1,1,1,1])
#h5.write_to_h5.VMD('segid',['s0','s0','s0','s0','s0','s1','s1','s1','s1','s1'])
#h5.write_to_h5.VMD('chain',['A','B'])
#h5.write_to_h5.VMD('atomicnumber',[0,1,2])

### Close file
h5.close()



#----------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------TIME INDEPENDENT----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#Create h5-file
h5_time_independent=h5md.h5md("File_time_independent.h5",system)

#Set arbitrary values for writing
for i in range(n_part):
    system.part[i].pos=np.array([0.1*(i+j),0.2*(i+j),0.3*(i+j)])
      
#Write to h5-file
h5_time_independent.write_to_h5.pos()                                         #With no additional parameters, use ()
h5_time_independent.write_to_h5.v(-1,"particles/atoms/velocity/","value")     #With additional parameters, use "-1" as first parameter 

#Close file
h5_time_independent.close()      
      
      