import h5py
import sys
import numpy as np


import espressomd._system as es
import espressomd
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate



include "myconfig.pxi"



class h5md(object):
  def __init__(self,filename,system):  
    self.system=system 
    self.filename=filename         
    self.h5_file=self.OpenFile(filename,"a")
    self.write_to_h5=self.write_to_h5(self)
    self.read_from_h5=self.read_from_h5(self)
    
  def OpenFile(self,filename,accesstype):
    file = h5py.File(filename,accesstype)
    return file
  
  def CreateDataset(self,file,dataset,shape,Maxshape,Dtype):
    dset = file.create_dataset(dataset,shape, maxshape=Maxshape, dtype=Dtype)
    return dset
  
  def ReadDataset(self,file,dataset_group,dataset_name):
    group=file[dataset_group]
    return group[dataset_name]
             
  def WriteValue(self,timestep,particle_id,value,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype,case,feature=1):
      if feature != 1:
        print "ERROR: Some necessary Espresso features for values used in h5-file are not activated"
        sys.exit()    
      try:
        self.dataset=self.CreateDataset(self.h5_file,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype)
      except:
        self.dataset=self.h5_file[h5_datasetpath]   
      if case=='step': 
        if(self.dataset.len()<=timestep+1):
          self.dataset.resize((timestep+1,1))
        self.dataset[timestep]=value  
      if case=='time': 
        if(self.dataset.len()<=timestep+1):
          self.dataset.resize((timestep+1,1))
        self.dataset[timestep]=value
        
      if case=='particle_time_independent': 
        n_part=self.system.n_part
        if(self.dataset.len()<=n_part+1):
          self.dataset.resize((n_part,self.dataset.shape[1]))
        self.dataset[particle_id]=value
      if case=='particle_time_dependent': 
        n_part=self.system.n_part
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,n_part,self.dataset.shape[2]))
        self.dataset[timestep,particle_id]=value 

      if case=='analyze_time_independent': 
        value_temp=[]
        try:  #value is array
          value_length=len(value)
          value_temp=value
        except:  #value is single
          value_length=1 
          value_temp.append(value)
        self.dataset.resize((value_length,))
        for i in range(0,value_length):
          self.dataset[i]=value_temp[i]
      if case=='analyze_time_dependent': 
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,self.dataset.shape[1]))
        self.dataset[timestep]=value 
               
      if case=='userdefined_value_time_independent':
        value_temp=[]
        try:  #value is array
          value_length=len(value)
          value_temp=value
        except:  #value is single
          value_length=1 
          value_temp.append(value)
        self.dataset.resize((value_length,))
        for i in range(0,value_length):
          self.dataset[i]=value_temp[i]
      if case=='userdefined_value_time_dependent': 
        try:  #value is array
          value_length=len(value)
        except: #value is single
          value_length=1 
        self.dataset.resize((timestep+1,value_length))
        self.dataset[timestep]=value

        
      if case=='box_edges_time_dependent': 
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,3,3))
        self.dataset[timestep,0,0]=self.system.box_l[0]
        self.dataset[timestep,1,1]=self.system.box_l[1]
        self.dataset[timestep,2,2]=self.system.box_l[2] 
      if case=='box_dimension_time_dependent':      
        self.dataset[0]=3
      if case=='box_boundary_time_dependent': 
        if(self.system.periodicity[0]==0):
          self.dataset[0]="none"
        else:
          self.dataset[0]="periodic"
        if(self.system.periodicity[1]==0):
          self.dataset[1]="none"
        else:
          self.dataset[1]="periodic"
        if(self.system.periodicity[2]==0):
          self.dataset[2]="none"
        else:
          self.dataset[2]="periodic" 



                  

  
  
#WRITE CLASS 
  class write_to_h5(object):
    def __init__(self,h5md):
      self.h5md=h5md          
      self.n_part=self.h5md.system.n_part   
    #time
    def time(self,timestep=-1,groupname="particles/atoms/Time/",datasetname="time"):
      self.h5md.WriteValue(timestep,-1,self.h5md.system.time,groupname+"/"+datasetname,(1,1),(None,1),'f8','time',1)
    #time step
    def time_step(self,timestep=-1,groupname="particles/atoms/Step/",datasetname="step"):
      self.h5md.WriteValue(timestep,-1,timestep,groupname+"/"+datasetname,(1,1),(None,1),'int64','step',1)

#PARTICLE PROPERTIES
    #position
    def pos(self,timestep=-1,groupname="particles/atoms/position/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',1)                      
    #image
    def image(self,timestep=-1,groupname="particles/atoms/image/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',1)             
    #velocity
    def v(self,timestep=-1,groupname="particles/atoms/velocity/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].v,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].v,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',1)      
    #force
    def f(self,timestep=-1,groupname="particles/atoms/force/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].f,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].f,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',1)       
    #bonds
    def bonds(self,groupname="particles/atoms/bond_from/XXXXXXXXX",datasetname="value"):
      self.n_part=self.h5md.system.n_part
      for i in range(0,self.n_part):
        for bond in self.h5md.system.part[i].bonds: 
          print("BOND:",bond[0],bond[1])         
          self.h5md.WriteValue(-1,i,i,"particles/"+groupname+"/bond_from/value",(1,1),(None,1),'int64','particle_time_independent',1) 
          self.h5md.WriteValue(-1,i,bond[1],"particles/"+groupname+"/bond_to/value",(1,1),(None,1),'int64','particle_time_independent',1) 
    #species
    def type(self,timestep=-1,groupname="particles/atoms/species/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].type,groupname+"/"+datasetname,(1,1),(None,1),'int64','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].type,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','particle_time_dependent',1)  
    #id
    def id(self,timestep=-1,groupname="particles/atoms/id/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].id,groupname+"/"+datasetname,(1,1),(None,1),'int64','particle_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].id,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','particle_time_dependent',1)
    #mass
    def mass(self,timestep=-1,groupname="particles/atoms/mass/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].mass,groupname+"/"+datasetname,(1,1),(None,1),'f8','particle_time_independent',MASS)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].mass,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','particle_time_dependent',MASS)        
    #omega_lab
    def omega_lab(self,timestep=-1,groupname="particles/atoms/omega_lab/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].omega_lab,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].omega_lab,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',ROTATION)        
    #rinertia
    def rinertia(self,timestep=-1,groupname="particles/atoms/rinertia/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].rinertia,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',ROTATIONAL_INERTIA)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].rinertia,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',ROTATIONAL_INERTIA)        
    #omega_body
    def omega_body(self,timestep=-1,groupname="particles/atoms/omega_body/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].omega_body,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].omega_body,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',ROTATION)        
    #torque_lab
    def torque_lab(self,timestep=-1,groupname="particles/atoms/torque_lab/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].torque_lab,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].torque_lab,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',ROTATION)        
    #quat
    def quat(self,timestep=-1,groupname="particles/atoms/quat/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].quat,groupname+"/"+datasetname,(1,4),(None,4),'f8','particle_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].quat,groupname+"/"+datasetname,(1,1,4),(None,None,4),'f8','particle_time_dependent',ROTATION)           
    #charge   
    def q(self,timestep=-1,groupname="particles/atoms/charge/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].q,groupname+"/"+datasetname,(1,1),(None,1),'f8','particle_time_independent',ELECTROSTATICS)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].q,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','particle_time_dependent',ELECTROSTATICS)               
    #virtual
    def virtual(self,timestep=-1,groupname="particles/atoms/virtual/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].virtual,groupname+"/"+datasetname,(1,1),(None,1),'int64','particle_time_independent',VIRTUAL_SITES_COM)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].virtual,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','particle_time_dependent',VIRTUAL_SITES_COM)     
    #vs_relative
    def vs_relative(self,timestep=-1,groupname="particles/atoms/vs_relative/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].vs_relative,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',VIRTUAL_SITES_RELATIVE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].vs_relative,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',VIRTUAL_SITES_RELATIVE) 
    #dipole
    def dip(self,timestep=-1,groupname="particles/atoms/dipole/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].dip,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',DIPOLES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].dip,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',DIPOLES)     
    #dipole_magnitude
    def dipm(self,timestep=-1,groupname="particles/atoms/dipole_magnitude/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].dipm,groupname+"/"+datasetname,(1,1),(None,1),'f8','particle_time_independent',DIPOLES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].dipm,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','particle_time_dependent',DIPOLES)     
    #external force
    def ext_force(self,timestep=-1,groupname="particles/atoms/ext_force/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].ext_force,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',EXTERNAL_FORCES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].ext_force,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',EXTERNAL_FORCES)    
    #external force particle fix
    def fix(self,timestep=-1,groupname="particles/atoms/ext_force_fix/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].fix,groupname+"/"+datasetname,(1,3),(None,3),'int64','particle_time_independent',EXTERNAL_FORCES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].fix,groupname+"/"+datasetname,(1,1,3),(None,None,3),'int64','particle_time_dependent',EXTERNAL_FORCES)   
    #external torque
    def ext_torque(self,timestep=-1,groupname="particles/atoms/ext_torque/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].ext_torque,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].ext_torque,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',ROTATION)        
    #gamma
    def gamma(self,timestep=-1,groupname="particles/atoms/gamma/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].gamma,groupname+"/"+datasetname,(1,1),(None,1),'f8','particle_time_independent',LANGEVIN_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].gamma,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','particle_time_dependent',LANGEVIN_PER_PARTICLE)         
    #temperature   
    def temp(self,timestep=-1,groupname="particles/atoms/temp/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].temp,groupname+"/"+datasetname,(1,1),(None,1),'f8','particle_time_independent',LANGEVIN_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].temp,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','particle_time_dependent',LANGEVIN_PER_PARTICLE)        
    #rotation   
    def rotation(self,timestep=-1,groupname="particles/atoms/rotation/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].rotation,groupname+"/"+datasetname,(1,1),(None,1),'int64','particle_time_independent',ROTATION_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].rotation,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','particle_time_dependent',ROTATION_PER_PARTICLE)                   
       
    
#OBSERVABLES
    #energy   
    def analyze_energy(self,timestep=-1,groupname="observables/energy/",datasetname="value"):
      if timestep == -1:
        for key in analyze.energy(self.h5md.system).keys():
          self.h5md.WriteValue(timestep,-1,analyze.energy(self.h5md.system)[key],groupname+"/"+datasetname+"_"+str(key),(1,),(None,),'f8','analyze_time_independent',1)          
      else:
        for key in analyze.energy(self.h5md.system).keys():
          self.h5md.WriteValue(timestep,-1,analyze.energy(self.h5md.system)[key],groupname+"/"+datasetname+"_"+str(key),(1,1),(None,None),'f8','analyze_time_dependent',1)  
    #analyze_linear_momentum   
    def analyze_linear_momentum(self,timestep=-1,groupname="observables/linear_momentum/",datasetname="value"):
      if timestep == -1:
        self.h5md.WriteValue(timestep,-1,analyze.analyze_linear_momentum(self.h5md.system),groupname+"/"+datasetname,(3,),(None,),'f8','analyze_time_independent',1)          
      else:
        self.h5md.WriteValue(timestep,-1,analyze.analyze_linear_momentum(self.h5md.system),groupname+"/"+datasetname,(1,3),(None,3),'f8','analyze_time_dependent',1)         
        


#BOX
    #box
    def box_edges(self,timestep=-1,groupname="particles/atoms/box/",datasetname="value"): 
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.box_l,groupname+"/"+datasetname,(1,1),(None,3),'f8','box_edges_time_independent',1)        
      else:
        self.h5md.WriteValue(timestep,-1,self.h5md.system.box_l,groupname+"/"+datasetname,(1,1,1),(None,3,3),'f8','box_edges_time_dependent')     
#     def box_boundary(self,timestep=-1,groupname="particles/atoms/box/",datasetname="boundary"): 
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,self.h5md.system.periodicity,groupname+"/"+datasetname,(1,3),(None,1),'f8','box_boundary_time_independent',1)        
#       else:     
#         self.h5md.WriteValue(timestep,-1,self.h5md.system.periodicity,groupname+"/"+datasetname,(1,3),(None,1),'S30','box_boundary_time_dependent')       
#     def box_dimension(self,timestep=-1,groupname="particles/atoms/box/",datasetname="dimension"):  
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,3,groupname+"/"+datasetname,(1,3),(None,3),'f8','box_dimension_time_independent',1)        
#       else:    
#         self.h5md.WriteValue(timestep,-1,3,groupname+"/"+datasetname,(1,1),(None,1),'int64','box_dimension_time_dependent') 
    #user defined dataset  
    def userdefined(self,timestep=-1,value="",groupname="",datasetname="",datatype='f8'):
      if timestep == -1:
        self.h5md.WriteValue(timestep,-1,value,groupname+"/"+datasetname,(1,),(None,),datatype,'userdefined_value_time_independent',1)        
      else:
        self.h5md.WriteValue(timestep,-1,value,groupname+"/"+datasetname,(1,1),(None,None),datatype,'userdefined_value_time_dependent',1) 


















   
      
      
      





























     
                      
#READ CLASS   
  class read_from_h5(object):
    def __init__(self,h5md):
      self.h5md=h5md

    #time
    def time(self,timestep,groupname="particles/atoms/Time/",datasetname="time"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      for i in range(self.h5md.value_dataset.shape[0]):
        self.h5md.system.time = self.h5md.value_dataset[i]  
    #time step
    def time_step(self,timestep,groupname="particles/atoms/Step/",datasetname="step"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      for i in range(self.h5md.value_dataset.shape[0]):
        timestep = self.h5md.value_dataset[i]  
    #Position
    def pos(self,timestep,groupname="particles/atoms/position/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      #Write positions in Espresso
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].pos = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].pos = self.h5md.value_dataset[timestep,i]                      
#     #image
#     def image(self,timestep=-1,groupname="particles/atoms/image/",datasetname="value"):
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,3),(None,3),'f8','particle_time_independent',1)        
#       else:
#         self.n_part=self.h5md.system.n_part
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','particle_time_dependent',1)             
    #velocity 
    def v(self,timestep,groupname="particles/atoms/velocity/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].v = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].v = self.h5md.value_dataset[timestep,i]    
    #force   
    def f(self,timestep,groupname="particles/atoms/force/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].f = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].f = self.h5md.value_dataset[timestep,i]   
#     #bonds
#     def bonds(self,groupname="particles/atoms/bond_from/XXXXXXXXX",datasetname="value"):
#       self.n_part=self.h5md.system.n_part
#       for i in range(0,self.n_part):
#         for bond in self.h5md.system.part[i].bonds: 
#           print("BOND:",bond[0],bond[1])         
#           self.h5md.WriteValue(-1,i,i,"particles/"+groupname+"/bond_from/value",(1,1),(None,1),'int64','particle_time_independent',1) 
#           self.h5md.WriteValue(-1,i,bond[1],"particles/"+groupname+"/bond_to/value",(1,1),(None,1),'int64','particle_time_independent',1) 
    #species  
    def type(self,timestep,groupname="particles/atoms/species/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].type = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].type = self.h5md.value_dataset[timestep,i][0]
    #id
    def id(self,timestep,groupname="particles/atoms/id/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2  :
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].id = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].id = self.h5md.value_dataset[timestep,i][0]  
    def mass(self,timestep,groupname="particles/atoms/mass/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].mass = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].mass = self.h5md.value_dataset[timestep,i][0]      
    #omega_lab 
    def omega_lab(self,timestep,groupname="particles/atoms/omega_lab/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].omega_lab = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].omega_lab = self.h5md.value_dataset[timestep,i]      
    #rinertia
    def rinertia(self,timestep,groupname="particles/atoms/rinertia/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].rinertia = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].rinertia = self.h5md.value_dataset[timestep,i]      
    #omega_body 
    def omega_body(self,timestep,groupname="particles/atoms/omega_body/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].omega_body = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].omega_body = self.h5md.value_dataset[timestep,i]      
    #torque_lab 
    def torque_lab(self,timestep,groupname="particles/atoms/torque_lab/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].torque_lab = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].torque_lab = self.h5md.value_dataset[timestep,i]      
    #quat     
    def quat(self,timestep,groupname="particles/atoms/quat/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].quat = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].quat = self.h5md.value_dataset[timestep,i]       
    #charge   
    def q(self,timestep,groupname="particles/atoms/charge/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].q = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].q = self.h5md.value_dataset[timestep,i][0]             
    #virtual
    def virtual(self,timestep,groupname="particles/atoms/virtual/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].virtual = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].virtual = self.h5md.value_dataset[timestep,i][0]     
    #vs_relative
    def vs_relative(self,timestep,groupname="particles/atoms/vs_relative/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].vs_relative = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].vs_relative = self.h5md.value_dataset[timestep,i] 
    #dipole 
    def dip(self,timestep,groupname="particles/atoms/dipole/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].dip = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].dip = self.h5md.value_dataset[timestep,i]   
    #dipole_magnitude    
    def dipm(self,timestep,groupname="particles/atoms/dipole_magnitude/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].dipm = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].dipm = self.h5md.value_dataset[timestep,i][0] 
    #external force
    def ext_force(self,timestep,groupname="particles/atoms/ext_force/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].ext_force = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].ext_force = self.h5md.value_dataset[timestep,i]    
    #external force particle fix
    def fix(self,timestep,groupname="particles/atoms/ext_force_fix/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].fix = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].fix = self.h5md.value_dataset[timestep,i]    
    #external torque
    def ext_torque(self,timestep,groupname="particles/atoms/ext_torque/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].ext_torque = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].ext_torque = self.h5md.value_dataset[timestep,i]       
    #gamma     
    def gamma(self,timestep,groupname="particles/atoms/gamma/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].gamma = self.h5md.value_dataset[i][0] 
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].gamma = self.h5md.value_dataset[timestep,i][0]
    #temperature    
    def temp(self,timestep,groupname="particles/atoms/temp/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].temp = self.h5md.value_dataset[i][0]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].temp = self.h5md.value_dataset[timestep,i][0]      
    #rotation      
    def rotation(self,timestep,groupname="particles/atoms/rotation/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          self.h5md.system.part[i].rotation = self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          self.h5md.system.part[i].rotation = self.h5md.value_dataset[timestep,i]                 
#     #box
#     def box_edges(self,timestep=-1,groupname="particles/atoms/box/",datasetname="value"): 
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,self.h5md.system.box_l,groupname+"/"+datasetname,(1,1),(None,3),'f8','box_edges_time_independent',1)        
#       else:
#         self.h5md.WriteValue(timestep,-1,self.h5md.system.box_l,groupname+"/"+datasetname,(1,1,1),(None,3,3),'f8','box_edges_time_dependent')  
#     def box_edges(self,timestep,groupname="particles/atoms/box/",datasetname="value"):
#       try:
#         self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
#       except:
#         print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
#         sys.exit()
#       if len(self.h5md.value_dataset.shape) == 2:
#         for i in range(self.h5md.value_dataset.shape[0]):
#           self.h5md.system.part[i].box_l = self.h5md.value_dataset[i]  
#       else:
#         for i in range(self.h5md.value_dataset.shape[1]):
#           self.h5md.system.part[i].box_l = self.h5md.value_dataset[timestep,i]    
#     def box_boundary(self,timestep=-1,groupname="particles/atoms/box/",datasetname="boundary"): 
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,self.h5md.system.periodicity,groupname+"/"+datasetname,(1,3),(None,1),'f8','box_boundary_time_independent',1)        
#       else:     
#         self.h5md.WriteValue(timestep,-1,self.h5md.system.periodicity,groupname+"/"+datasetname,(1,3),(None,1),'S30','box_boundary_time_dependent') 
#     def XXX(self,timestep,groupname="particles/atoms/XXX/",datasetname="value"):
#       try:
#         self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
#       except:
#         print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
#         sys.exit()
#       if len(self.h5md.value_dataset.shape) == 2:
#         for i in range(self.h5md.value_dataset.shape[0]):
#           self.h5md.system.part[i].XXX = self.h5md.value_dataset[i]  
#       else:
#         for i in range(self.h5md.value_dataset.shape[1]):
#           self.h5md.system.part[i].XXX = self.h5md.value_dataset[timestep,i]      
#     def box_dimension(self,timestep=-1,groupname="particles/atoms/box/",datasetname="dimension"):  
#       if timestep == -1:
#         for i in range(0,self.n_part):
#           self.h5md.WriteValue(-1,i,3,groupname+"/"+datasetname,(1,3),(None,3),'f8','box_dimension_time_independent',1)        
#       else:    
#         self.h5md.WriteValue(timestep,-1,3,groupname+"/"+datasetname,(1,1),(None,1),'int64','box_dimension_time_dependent') 
#     def XXX(self,timestep,groupname="particles/atoms/XXX/",datasetname="value"):
#       try:
#         self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
#       except:
#         print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
#         sys.exit()
#       if len(self.h5md.value_dataset.shape) == 2:
#         for i in range(self.h5md.value_dataset.shape[0]):
#           self.h5md.system.part[i].XXX = self.h5md.value_dataset[i]  
#       else:
#         for i in range(self.h5md.value_dataset.shape[1]):
#           self.h5md.system.part[i].XXX = self.h5md.value_dataset[timestep,i] 
    #user defined dataset   
    def userdefined(self,timestep,groupname="particles/atoms/userdefined/",datasetname="value"):
      try:
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,groupname,datasetname)
      except:
        print "Error: No "+groupname+"/"+datasetname+" dataset in h5-file exisitng"
        sys.exit()
      if len(self.h5md.value_dataset.shape) == 2:
        for i in range(self.h5md.value_dataset.shape[0]):
          return self.h5md.value_dataset[i]  
      else:
        for i in range(self.h5md.value_dataset.shape[1]):
          return self.h5md.value_dataset[timestep,i] 