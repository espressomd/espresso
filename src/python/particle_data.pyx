#from espresso.utils cimport ERROR 
cimport numpy as np
import numpy as np
cimport utils


cdef class ParticleHandle:
  def __cinit__(self, _id):
#    utils.init_intlist(self.particleData.el)
    utils.init_intlist(&(self.particleData.bl))
    self.id=_id

  cdef int update_particle_data(self) except -1:
#    utils.realloc_intlist(self.particleData.el, 0)
    utils.realloc_intlist(&(self.particleData.bl), 0)
      
    if get_particle_data(self.id, &self.particleData):
      raise Exception("Error updating particle data")
    else: 
      return 0

  property type:
    def __set__(self, _type):
      if isinstance(_type, int) and _type >= 0:  
        if set_particle_type(self.id, _type) == 1:
          print 'set particle position first'
      else:
        print 'type must be an integer >= 0'
    def __get__(self):
      self.update_particle_data()
      return self.particleData.p.type

  property pos:
    def __set__(self, _pos):
      cdef double mypos[3]
      for i in range(3):
        if not isinstance(_pos[i], float):
          print 'position must be float'
        else:
          mypos[i]=_pos[i]
      if place_particle(self.id, mypos) == -1:
        print 'particle could not be set'
    def __get__(self):
      self.update_particle_data()
      return np.array([self.particleData.r.p[0],\
                       self.particleData.r.p[1],\
                       self.particleData.r.p[2]])

  property v:
    def __set__(self, _v):
      cdef double myv[3]
      for i in range(3):
        if not isinstance(_v[i], float):
          print 'velocity must be float'
        else:
          myv[i]=_v[i]
      if set_particle_v(self.id, myv) == 1:
        print 'set particle position first'
    def __get__(self):
      self.update_particle_data()
      return np.array([ self.particleData.m.v[0],\
                        self.particleData.m.v[1],\
                        self.particleData.m.v[2]])

  property f:
    def __set__(self, _f):
      cdef double myf[3]
      for i in range(3):
        if not isinstance(_f[i], float):
          print 'force must be float'
        else:
          myf[i]=_f[i]
      if set_particle_f(self.id, myf) == 1:
        print 'set particle position first'
    def __get__(self):
      self.update_particle_data()
      return np.array([ self.particleData.f.f[0],\
                        self.particleData.f.f[1],\
                        self.particleData.f.f[2]])


cdef class particleList:
  def __getitem__(self, key):
    return ParticleHandle(key)


