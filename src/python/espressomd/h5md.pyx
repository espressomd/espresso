import h5py
import sys
import numpy as np
include "myconfig.pxi"

class h5md(object):
  def __init__(self,filename,system):  
    self.system=system 
    self.filename=filename         
    self.h5_file=self.OpenFile(filename,"a")
    self.h5_write_particles=self.h5_write_particles(self)
    self.h5_read_particles=self.h5_read_particles(self)
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
             
  def WriteValue(self,timestep,particle_id,particle_value,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype,case,feature=1):
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
        self.dataset[timestep]=timestep  
      if case=='time': 
        if(self.dataset.len()<=timestep+1):
          self.dataset.resize((timestep+1,1))
        self.dataset[timestep]=self.system.time
      if case=='value': 
        n_part=self.system.n_part
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,n_part,self.dataset.shape[2]))
        self.dataset[timestep,particle_id]=particle_value 
      if case=='value_scalar': 
        n_part=self.system.n_part
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,n_part,1))
        self.dataset[timestep,particle_id]=particle_value 
      if case=='value_fix': 
        n_part=self.system.n_part
        if(self.dataset.len()<=n_part+1):
          self.dataset.resize((n_part,self.dataset.shape[1]))
        self.dataset[particle_id]=particle_value
      if case=='box_dimension':      
        self.dataset[0]=3
      if case=='box_periodicity': 
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
      if case=='box_edges': 
        if(self.dataset.len()<=timestep+1): 
          self.dataset.resize((timestep+1,3,3))
        #TODO
        self.dataset[timestep,0,0]=self.system.box_l[0]
        self.dataset[timestep,1,1]=self.system.box_l[1]
        self.dataset[timestep,2,2]=self.system.box_l[2]  
        
  def WriteObservableValue(self,timestep,value,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype,case):   
    try:
      self.dataset=self.CreateDataset(self.h5_file,h5_datasetpath,h5_shape,h5_maxshape,h5_Dtype)
    except:
      self.dataset=self.h5_file[h5_datasetpath]   
    if case=='step': 
      if(self.dataset.len()<=timestep+1):
        self.dataset.resize((timestep+1,1))
      self.dataset[timestep]=timestep  
    if case=='time': 
      if(self.dataset.len()<=timestep+1):
        self.dataset.resize((timestep+1,1))
      self.dataset[timestep]=self.system.time
    if case=='value': 
      try:
        value_length=len(value)
      except:
        value_length=1 
      self.dataset.resize((timestep+1,value_length))
      for i in range(0,value_length):
        self.dataset[timestep,i]=value[i]
  
#WRITE CLASS 
  class h5_write_particles(object):
    def __init__(self,self_h5md_class):
      self.self_h5md_class=self_h5md_class          
      self.n_part=self.self_h5md_class.system.n_part   
    #position
    def pos(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/position/step",(1,1),(None,1),'int64','step')
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/position/time",(1,1),(None,1),'f8','time')
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].pos,"particles/"+groupname+"/position/value",(1,1,3),(None,None,3),'f8','value')             
    #pmage
    def image(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/image/step",(1,1),(None,1),'int64','step')
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/image/time",(1,1),(None,1),'f8','time')
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].pos,"particles/"+groupname+"/image/value",(1,1,3),(None,None,3),'f8','value')             
    #velocity
    def v(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/velocity/step",(1,1),(None,1),'int64','step')
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/velocity/time",(1,1),(None,1),'f8','time')
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].v,"particles/"+groupname+"/velocity/value",(1,1,3),(None,None,3),'f8','value')      
    #force
    def f(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/force/step",(1,1),(None,1),'int64','step')
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/force/time",(1,1),(None,1),'f8','time')
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].f,"particles/"+groupname+"/force/value",(1,1,3),(None,None,3),'f8','value')       
    #species
    def type(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].type,"particles/"+groupname+"/species/value",(1,1),(None,1),'int64','value_fix')  
    #id
    def id(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].id,"particles/"+groupname+"/id/value",(1,1),(None,1),'int64','value_fix')
    #mass
    def mass(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].mass,"particles/"+groupname+"/mass/value",(1,1),(None,1),'f8','value_fix',MASS)
    #omega_lab
    def omega_lab(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/omega_lab/step",(1,1),(None,1),'int64','step',ROTATION)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/omega_lab/time",(1,1),(None,1),'f8','time',ROTATION)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].omega_lab,"particles/"+groupname+"/omega_lab/value",(1,1,3),(None,None,3),'f8','value',ROTATION)        
    #rinertia
    def rinertia(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/rinertia/step",(1,1),(None,1),'int64','step',ROTATIONAL_INERTIA)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/rinertia/time",(1,1),(None,1),'f8','time',ROTATIONAL_INERTIA)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].rinertia,"particles/"+groupname+"/rinertia/value",(1,1,3),(None,None,3),'f8','value',ROTATIONAL_INERTIA)        
    #omega_body
    def omega_body(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/omega_body/step",(1,1),(None,1),'int64','step',ROTATION)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/omega_body/time",(1,1),(None,1),'f8','time',ROTATION)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].omega_body,"particles/"+groupname+"/omega_body/value",(1,1,3),(None,None,3),'f8','value',ROTATION)        
    #torque_lab
    def torque_lab(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/torque_lab/step",(1,1),(None,1),'int64','step',ROTATION)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/torque_lab/time",(1,1),(None,1),'f8','time',ROTATION)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].torque_lab,"particles/"+groupname+"/torque_lab/value",(1,1,3),(None,None,3),'f8','value',ROTATION)        
    #quat
    def quat(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/quat/step",(1,1),(None,1),'int64','step',ROTATION)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/quat/time",(1,1),(None,1),'f8','time',ROTATION)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].quat,"particles/"+groupname+"/quat/value",(1,1,4),(None,None,4),'f8','value',ROTATION)            
    #charge
    def q(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].q,"particles/"+groupname+"/charge/value",(1,1),(None,1),'int64','value_fix',ELECTROSTATICS)   
    #virtual
    def virtual(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].virtual,"particles/"+groupname+"/virtual/value",(1,1),(None,1),'int64','value_fix',VIRTUAL_SITES_COM)
    #vs_relative
    def vs_relative(self,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(-1,i,self.self_h5md_class.system.part[i].vs_relative,"particles/"+groupname+"/vs_relative/value",(1,2),(None,2),'f8','value_fix',VIRTUAL_SITES_RELATIVE)
    #cipole
    def dip(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/dipole/step",(1,1),(None,1),'int64','step',DIPOLES)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/dipole/time",(1,1),(None,1),'f8','time',DIPOLES)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].dip,"particles/"+groupname+"/dipole/value",(1,1,3),(None,None,3),'f8','value',DIPOLES)     
    #dipole_magnitude
    def dipm(self,timestep,groupname="atoms"):
      self.n_part=self.self_h5md_class.system.n_part
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/dipole_magnitude/step",(1,1),(None,1),'int64','step',DIPOLES)
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/dipole_magnitude/time",(1,1),(None,1),'f8','time',DIPOLES)
      for i in range(0,self.n_part):
        self.self_h5md_class.WriteValue(timestep,i,self.self_h5md_class.system.part[i].dipm,"particles/"+groupname+"/dipole_magnitude/value",(1,1,3),(None,None,3),'f8','value',DIPOLES)     
    #box
    def box(self,timestep=0,groupname="atoms"): 
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/box/step",(1,1),(None,1),'int64','step')
      self.self_h5md_class.WriteValue(timestep,-1,-1,"particles/"+groupname+"/box/time",(1,1),(None,1),'f8','time')
      self.self_h5md_class.WriteValue(timestep,-1,self.self_h5md_class.system.box_l,"particles/"+groupname+"/box/value",(1,1,3),(None,None,3),'f8','box_edges')     
      self.self_h5md_class.WriteValue(timestep,-1,self.self_h5md_class.system.periodicity,"particles/"+groupname+"/box/boundary",(3,1),(3,1),'S30','box_periodicity')       
      self.self_h5md_class.WriteValue(timestep,-1,3,"particles/"+groupname+"/box/dimension",(1,1),(None,1),'int64','box_dimension') 
                      
#READ CLASS   
  class h5_read_particles(object):
    def __init__(self,self_h5md_class):
      self.self_h5md_class=self_h5md_class
     
    #Position
    def pos(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/position','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/position','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/position','value')
      except:
        print "Error: No particles/"+groupname+"/position/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].pos = self.self_h5md_class.value_dataset[timestep,i]                      
    #Image
    def image(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/image','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/image','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/image','value')
      except:
        print "Error: No particles/"+groupname+"/image/value,time,step dataset in h5-file available"
        sys.exit()
      #Write image and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].pos = self.self_h5md_class.value_dataset[timestep,i]                
    #Velocity
    def v(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/velocity','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/velocity','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/velocity','value')
      except:
        print "Error: No particles/"+groupname+"/velocity/value,time,step dataset in h5-file available"
        sys.exit()
      #Write velocity and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].v = self.self_h5md_class.value_dataset[timestep,i]                
    #Force
    def f(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/force','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/force','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/force','value')
      except:
        print "Error: No particles/"+groupname+"/force/value,time,step dataset in h5-file available"
        sys.exit()
      #Write force and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].f = self.self_h5md_class.value_dataset[timestep,i]               
    #Species
    def type(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/species','value')
      except:
        print "Error: No particles/"+groupname+"/species/value dataset in h5-file available"
        sys.exit()
      #Write type in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].type = int(self.self_h5md_class.dataset[i])                
    #ID
    def id(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/id','value')
      except:
        print "Error: No particles/"+groupname+"/id/value dataset in h5-file available"
        sys.exit()
      #Write ID in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].id = int(self.self_h5md_class.dataset[i])             
    #Mass
    def mass(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/mass','value')
      except:
        print "Error: No particles/"+groupname+"/mass/value dataset in h5-file available"
        sys.exit()
      #Write mass in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].mass = int(self.self_h5md_class.dataset[i])     
    #omega_lab
    def omega_lab(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_lab','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_lab','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_lab','value')
      except:
        print "Error: No particles/"+groupname+"/omega_lab/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].omega_lab = self.self_h5md_class.value_dataset[timestep,i]  
    #rinertia
    def rinertia(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/rinertia','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/rinertia','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/rinertia','value')
      except:
        print "Error: No particles/"+groupname+"/rinertia/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].rinertia = self.self_h5md_class.value_dataset[timestep,i]  
    #omega_body
    def omega_body(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_body','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_body','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/omega_body','value')
      except:
        print "Error: No particles/"+groupname+"/omega_body/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].omega_body = self.self_h5md_class.value_dataset[timestep,i]    
    #torque_lab
    def torque_lab(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/torque_lab','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/torque_lab','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/torque_lab','value')
      except:
        print "Error: No particles/"+groupname+"/torque_lab/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].torque_lab = self.self_h5md_class.value_dataset[timestep,i]     
    #quat
    def quat(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/quat','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/quat','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/quat','value')
      except:
        print "Error: No particles/"+groupname+"/quat/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].quat = self.self_h5md_class.value_dataset[timestep,i]         
    #charge
    def charge(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/charge','value')
      except:
        print "Error: No particles/"+groupname+"/charge/value dataset in h5-file available"
        sys.exit()
      #Write charge in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].charge = int(self.self_h5md_class.dataset[i])        
    #virtual
    def virtual(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/virtual','value')
      except:
        print "Error: No particles/"+groupname+"/virtual/value dataset in h5-file available"
        sys.exit()
      #Write virtual in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].virtual = int(self.self_h5md_class.dataset[i])
    #vs_relative
    def vs_relative(self,groupname="atoms"):
      try:
        self.self_h5md_class.dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/vs_relative','value')
      except:
        print "Error: No particles/"+groupname+"/vs_relative/value dataset in h5-file available"
        sys.exit()
      #Write vs_relative in Espresso
      for i in range(self.self_h5md_class.dataset.shape[0]):
        self.self_h5md_class.system.part[i].vs_relative = int(self.self_h5md_class.dataset[i])        
    #dipole
    def dipole(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipole','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipole','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipole','value')
      except:
        print "Error: No particles/"+groupname+"/dipole/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].dipole = self.self_h5md_class.value_dataset[timestep,i]            
    #dipoleMagnitude
    def dipoleMagnitude(self,timestep,groupname="atoms"):
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipoleMagnitude','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipoleMagnitude','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/dipoleMagnitude','value')
      except:
        print "Error: No particles/"+groupname+"/dipoleMagnitude/value,time,step dataset in h5-file available"
        sys.exit()
      #Write positions and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      for i in range(self.self_h5md_class.value_dataset.shape[1]):
        self.self_h5md_class.system.part[i].dipoleMagnitude = self.self_h5md_class.value_dataset[timestep,i]                           
    #Box
    def box(self,timestep=0,groupname="atoms"):
      try:
        self.self_h5md_class.dimension_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/box','dimension')
        self.self_h5md_class.boundary_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/box','boundary')
      except:
        print "Error: No particles/"+groupname+"/box/boundary,dimension dataset in h5-file available"
        sys.exit()
        
      if(self.self_h5md_class.boundary_dataset[0]=="none"):
        self.self_h5md_class.system.periodicity[0]=0
      else:
        self.self_h5md_class.system.periodicity[0]=1
      if(self.self_h5md_class.boundary_dataset[0]=="none"):
        self.self_h5md_class.system.periodicity[1]=0
      else:
        self.self_h5md_class.system.periodicity[1]=1
      if(self.self_h5md_class.boundary_dataset[0]=="none"):
        self.self_h5md_class.system.periodicity[2]=0
      else:
        self.self_h5md_class.system.periodicity[2]=1
             
      try:
        self.self_h5md_class.step_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/box/edges','step')
        self.self_h5md_class.time_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/box/edges','time')
        self.self_h5md_class.value_dataset=self.self_h5md_class.ReadDataset(self.self_h5md_class.h5_file,'particles/'+groupname+'/box/edges','value')
      except:
        print "Error: No particles/"+groupname+"/box/edges/value,time,step dataset in h5-file available"
        sys.exit()
        
      #Write velocity and time in Espresso
      self.self_h5md_class.system.time = self.self_h5md_class.time_dataset[timestep]
      #TODO
      self.self_h5md_class.system.box_l[0] = self.self_h5md_class.value_dataset[timestep,0,0]      
      self.self_h5md_class.system.box_l[1] = self.self_h5md_class.value_dataset[timestep,1,1]      
      self.self_h5md_class.system.box_l[2] = self.self_h5md_class.value_dataset[timestep,2,2]  


#OBSERVABLES
  def h5_write_observable(self,timestep,value,observablename,groupname="observable"):
    self.WriteObservableValue(timestep,-1,"observables/"+groupname+"/"+observablename+"/step",(1,1),(None,1),'int64','step')
    self.WriteObservableValue(timestep,-1,"observables/"+groupname+"/"+observablename+"/time",(1,1),(None,1),'f8','time')
    self.WriteObservableValue(timestep,value,"observables/"+groupname+"/"+observablename+"/value",(1,1),(None,None),'f8','value')                      
  def h5_read_observable(self,timestep,observablename,groupname):  
    try:
      self.step_dataset=self.ReadDataset(self.h5_file,'observables/'+groupname+'/'+observablename,'step')
      self.time_dataset=self.ReadDataset(self.h5_file,'observables/'+groupname+'/'+observablename,'time')
      self.value_dataset=self.ReadDataset(self.h5_file,'observables/'+groupname+'/'+observablename,'value')
      return self.step_dataset,self.time_dataset,self.value_dataset
    except:
      print "Error: No observables/"+groupname+"/"+observablename+"/value,time,step dataset in h5-file available"
      sys.exit()     
  
  


