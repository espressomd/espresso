import h5py
import sys
import numpy as np
include "myconfig.pxi"

class h5md(object):
  def __init__(self,filename,system):  
    self.system=system 
    self.filename=filename         
    self.h5_file=self.OpenFile(filename,"a")
    self.write_to_h5=self.write_to_h5(self)
    self.read_from_h5=self.read_from_h5(self)
    #self.h5_write_vmd_parameters_extra=self.h5_write_vmd_parameters_extra(self)
    #self.h5_read_vmd_parameters_extra=self.h5_read_vmd_parameters_extra(self)
    
  def OpenFile(self,filename,accesstype):
    file = h5py.File(filename,accesstype)
    return file
  
  def DatasetSize(self,dataset):
    return self.file[dataset].shape
  
  def WriteAttributes(self,dataset,name,value):
    self.file[dataset].attrs[name] = value
  
  def ReadAttributes(self,dataset,name):
    return self.file[dataset].attrs[name]
  
  def CreateDataset(self,file,dataset,shape,Maxshape,Dtype):
    dset = file.create_dataset(dataset,shape, maxshape=Maxshape, dtype=Dtype)
    return dset
  
  def ReadDataset(self,file,dataset_group,dataset_name):
    group=file[dataset_group]
    return group[dataset_name]
             
  def WriteValue(self,timestep,particle_id,value_time_independent,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype,case,feature=1):
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
        self.dataset[timestep]=value_time_independent  
      if case=='time': 
        if(self.dataset.len()<=timestep+1):
          self.dataset.resize((timestep+1,1))
        self.dataset[timestep]=value_time_independent
      if case=='value_time_dependent': 
        n_part=self.system.n_part
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,n_part,self.dataset.shape[2]))
        self.dataset[timestep,particle_id]=value_time_independent 
      if case=='value_scalar': 
        n_part=self.system.n_part
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,n_part,1))
        self.dataset[timestep,particle_id]=value_time_independent 
      if case=='value_time_independent': 
        n_part=self.system.n_part
        if(self.dataset.len()<=n_part+1):
          self.dataset.resize((n_part,self.dataset.shape[1]))
        self.dataset[particle_id]=value_time_independent
        
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
          
      if case=='value_userdefined': 
        try:
          value_length=len(value_time_independent)
        except:
          value_length=1 
        self.dataset.resize((timestep+1,value_length))
        for i in range(0,value_length):
          self.dataset[timestep,i]=value_time_independent[i]
                
  
  
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
    #position
    def pos(self,timestep=-1,groupname="particles/atoms/position/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',1)                      
    #image
    def image(self,timestep=-1,groupname="particles/atoms/image/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].pos,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',1)             
    #velocity
    def v(self,timestep=-1,groupname="particles/atoms/velocity/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].v,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].v,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',1)      
    #force
    def f(self,timestep=-1,groupname="particles/atoms/force/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].f,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].f,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',1)       
    #bonds
    def bonds(self,groupname="particles/atoms/bond_from/XXXXXXXXX",datasetname="value"):
      self.n_part=self.h5md.system.n_part
      for i in range(0,self.n_part):
        for bond in self.h5md.system.part[i].bonds: 
          print("BOND:",bond[0],bond[1])         
          self.h5md.WriteValue(-1,i,i,"particles/"+groupname+"/bond_from/value",(1,1),(None,1),'int64','value_time_independent',1) 
          self.h5md.WriteValue(-1,i,bond[1],"particles/"+groupname+"/bond_to/value",(1,1),(None,1),'int64','value_time_independent',1) 
    #species
    def type(self,timestep=-1,groupname="particles/atoms/species/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].type,groupname+"/"+datasetname,(1,1),(None,1),'int64','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].type,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','value_time_dependent',1)  
    #id
    def id(self,timestep=-1,groupname="particles/atoms/id/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].id,groupname+"/"+datasetname,(1,1),(None,1),'int64','value_time_independent',1)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].id,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','value_time_dependent',1)
    #mass
    def mass(self,timestep=-1,groupname="particles/atoms/mass/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].mass,groupname+"/"+datasetname,(1,1),(None,1),'f8','value_time_independent',MASS)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].mass,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','value_time_dependent',MASS)        
    #omega_lab
    def omega_lab(self,timestep=-1,groupname="particles/atoms/omega_lab/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].omega_lab,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].omega_lab,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',ROTATION)        
    #rinertia
    def rinertia(self,timestep=-1,groupname="particles/atoms/rinertia/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].rinertia,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',ROTATIONAL_INERTIA)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].rinertia,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',ROTATIONAL_INERTIA)        
    #omega_body
    def omega_body(self,timestep=-1,groupname="particles/atoms/omega_body/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].omega_body,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].omega_body,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',ROTATION)        
    #torque_lab
    def torque_lab(self,timestep=-1,groupname="particles/atoms/torque_lab/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].torque_lab,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].torque_lab,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',ROTATION)        
    #quat
    def quat(self,timestep=-1,groupname="particles/atoms/quat/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].quat,groupname+"/"+datasetname,(1,4),(None,4),'f8','value_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].quat,groupname+"/"+datasetname,(1,1,4),(None,None,4),'f8','value_time_dependent',ROTATION)           
    #charge   
    def q(self,timestep=-1,groupname="particles/atoms/charge/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].q,groupname+"/"+datasetname,(1,1),(None,1),'f8','value_time_independent',ELECTROSTATICS)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].q,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','value_time_dependent',ELECTROSTATICS)               
    #virtual
    def virtual(self,timestep=-1,groupname="particles/atoms/virtual/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].virtual,groupname+"/"+datasetname,(1,1),(None,1),'int64','value_time_independent',VIRTUAL_SITES_COM)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].virtual,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','value_time_dependent',VIRTUAL_SITES_COM)     
    #vs_relative
    def vs_relative(self,timestep=-1,groupname="particles/atoms/vs_relative/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].vs_relative,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',VIRTUAL_SITES_RELATIVE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].vs_relative,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',VIRTUAL_SITES_RELATIVE) 
    #dipole
    def dip(self,timestep=-1,groupname="particles/atoms/dipole/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].dip,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',DIPOLES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].dip,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',DIPOLES)     
    #dipole_magnitude
    def dipm(self,timestep=-1,groupname="particles/atoms/dipole_magnitude/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].dipm,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',DIPOLES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].dipm,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',DIPOLES)     
    #external force
    def ext_force(self,timestep=-1,groupname="particles/atoms/ext_force/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].ext_force,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',EXTERNAL_FORCES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].ext_force,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',EXTERNAL_FORCES)    
    #external force particle fix
    def fix(self,timestep=-1,groupname="particles/atoms/ext_force_fix/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].fix,groupname+"/"+datasetname,(1,3),(None,3),'int64','value_time_independent',EXTERNAL_FORCES)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].fix,groupname+"/"+datasetname,(1,1,3),(None,None,3),'int64','value_time_dependent',EXTERNAL_FORCES)   
    #external torque
    def ext_torque(self,timestep=-1,groupname="particles/atoms/ext_torque/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].ext_torque,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',ROTATION)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].ext_torque,groupname+"/"+datasetname,(1,1,3),(None,None,3),'f8','value_time_dependent',ROTATION)        
    #gamma
    def gamma(self,timestep=-1,groupname="particles/atoms/gamma/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].gamma,groupname+"/"+datasetname,(1,1),(None,1),'f8','value_time_independent',LANGEVIN_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].gamma,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','value_time_dependent',LANGEVIN_PER_PARTICLE)         
    #temperature   
    def temp(self,timestep=-1,groupname="particles/atoms/temp/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].temp,groupname+"/"+datasetname,(1,1),(None,1),'f8','value_time_independent',LANGEVIN_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].temp,groupname+"/"+datasetname,(1,1,1),(None,None,1),'f8','value_time_dependent',LANGEVIN_PER_PARTICLE)        
    #rotation   
    def rotation(self,timestep=-1,groupname="particles/atoms/rotation/",datasetname="value"):
      if timestep == -1:
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].rotation,groupname+"/"+datasetname,(1,1),(None,1),'int64','value_time_independent',ROTATION_PER_PARTICLE)        
      else:
        self.n_part=self.h5md.system.n_part
        for i in range(0,self.n_part):
          self.h5md.WriteValue(timestep,i,self.h5md.system.part[i].rotation,groupname+"/"+datasetname,(1,1,1),(None,None,1),'int64','value_time_dependent',ROTATION_PER_PARTICLE)                   
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
        for i in range(0,self.n_part):
          self.h5md.WriteValue(-1,i,self.h5md.system.part[i].XXX,groupname+"/"+datasetname,(1,3),(None,3),'f8','value_time_independent',1)        
      else:
        self.h5md.WriteValue(timestep,-1,value,groupname+"/"+datasetname,(1,1),(None,None),datatype,'value_userdefined',1)       
 
      
      
      
      
      
                      
#READ CLASS   
  class read_from_h5(object):
    def __init__(self,h5md):
      self.h5md=h5md
     
    #Position
    def pos(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/position','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/position','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/position','value')
      except:
        print "Error: No particles/"+groupname+"/position/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].pos = self.h5md.value_dataset[timestep,i]                      
    #Image
    def image(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/image','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/image','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/image','value')
      except:
        print "Error: No particles/"+groupname+"/image/value,time,step dataset in h5-file available"
        sys.exit()
      #Write image and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].pos = self.h5md.value_dataset[timestep,i]                
    #Velocity
    def v(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/velocity','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/velocity','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/velocity','value')
      except:
        print "Error: No particles/"+groupname+"/velocity/value,time,step dataset in h5-file available"
        sys.exit()
      #Write velocity and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].v = self.h5md.value_dataset[timestep,i]                
    #Force
    def f(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/force','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/force','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/force','value')
      except:
        print "Error: No particles/"+groupname+"/force/value,time,step dataset in h5-file available"
        sys.exit()
      #Write force and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].f = self.h5md.value_dataset[timestep,i]               
    #Species
    def type(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/species','value')
      except:
        print "Error: No particles/"+groupname+"/species/value dataset in h5-file available"
        sys.exit()
      #Write type in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].type = int(self.h5md.dataset[i])                
    #ID
    def id(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/id','value')
      except:
        print "Error: No particles/"+groupname+"/id/value dataset in h5-file available"
        sys.exit()
      #Write ID in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].id = int(self.h5md.dataset[i])             
    #Mass
    def mass(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/mass','value')
      except:
        print "Error: No particles/"+groupname+"/mass/value dataset in h5-file available"
        sys.exit()
      #Write mass in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].mass = int(self.h5md.dataset[i])     
    #omega_lab
    def omega_lab(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_lab','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_lab','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_lab','value')
      except:
        print "Error: No particles/"+groupname+"/omega_lab/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].omega_lab = self.h5md.value_dataset[timestep,i]  
    #rinertia
    def rinertia(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/rinertia','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/rinertia','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/rinertia','value')
      except:
        print "Error: No particles/"+groupname+"/rinertia/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].rinertia = self.h5md.value_dataset[timestep,i]  
    #omega_body
    def omega_body(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_body','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_body','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/omega_body','value')
      except:
        print "Error: No particles/"+groupname+"/omega_body/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].omega_body = self.h5md.value_dataset[timestep,i]    
    #torque_lab
    def torque_lab(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/torque_lab','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/torque_lab','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/torque_lab','value')
      except:
        print "Error: No particles/"+groupname+"/torque_lab/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].torque_lab = self.h5md.value_dataset[timestep,i]     
    #quat
    def quat(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/quat','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/quat','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/quat','value')
      except:
        print "Error: No particles/"+groupname+"/quat/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].quat = self.h5md.value_dataset[timestep,i]         
    #charge
    def charge(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/charge','value')
      except:
        print "Error: No particles/"+groupname+"/charge/value dataset in h5-file available"
        sys.exit()
      #Write charge in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].charge = int(self.h5md.dataset[i])        
    #virtual
    def virtual(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/virtual','value')
      except:
        print "Error: No particles/"+groupname+"/virtual/value dataset in h5-file available"
        sys.exit()
      #Write virtual in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].virtual = int(self.h5md.dataset[i])
    #vs_relative
    def vs_relative(self,groupname="atoms"):
      try:
        self.h5md.dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/vs_relative','value')
      except:
        print "Error: No particles/"+groupname+"/vs_relative/value dataset in h5-file available"
        sys.exit()
      #Write vs_relative in Espresso
      for i in range(self.h5md.dataset.shape[0]):
        self.h5md.system.part[i].vs_relative = int(self.h5md.dataset[i])        
    #dipole
    def dipole(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipole','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipole','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipole','value')
      except:
        print "Error: No particles/"+groupname+"/dipole/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].dipole = self.h5md.value_dataset[timestep,i]            
    #dipoleMagnitude
    def dipoleMagnitude(self,timestep,groupname="atoms"):
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipoleMagnitude','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipoleMagnitude','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/dipoleMagnitude','value')
      except:
        print "Error: No particles/"+groupname+"/dipoleMagnitude/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      for i in range(self.h5md.value_dataset.shape[1]):
        self.h5md.system.part[i].dipoleMagnitude = self.h5md.value_dataset[timestep,i]                           
    #Box
    def box(self,timestep,groupname="atoms"):
      try:
        self.h5md.dimension_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/box','dimension')
        self.h5md.boundary_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/box','boundary')
      except:
        print "Error: No particles/"+groupname+"/box/boundary,dimension dataset in h5-file available"
        sys.exit()
        
      if(self.h5md.boundary_dataset[0]=="none"):
        self.h5md.system.periodicity[0]=0
      else:
        self.h5md.system.periodicity[0]=1
      if(self.h5md.boundary_dataset[0]=="none"):
        self.h5md.system.periodicity[1]=0
      else:
        self.h5md.system.periodicity[1]=1
      if(self.h5md.boundary_dataset[0]=="none"):
        self.h5md.system.periodicity[2]=0
      else:
        self.h5md.system.periodicity[2]=1
             
      try:
        self.h5md.step_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/box/edges','step')
        self.h5md.time_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/box/edges','time')
        self.h5md.value_dataset=self.h5md.ReadDataset(self.h5md.h5_file,'particles/'+groupname+'/box/edges','value')
      except:
        print "Error: No particles/"+groupname+"/box/edges/value,time,step dataset in h5-file available"
        sys.exit()
        
      #Write velocity and time in Espresso
      self.h5md.system.time = self.h5md.time_dataset[timestep]
      #TODO
      self.h5md.system.box_l[0] = self.h5md.value_dataset[timestep,0,0]      
      self.h5md.system.box_l[1] = self.h5md.value_dataset[timestep,1,1]      
      self.h5md.system.box_l[2] = self.h5md.value_dataset[timestep,2,2]  


    
  
  


