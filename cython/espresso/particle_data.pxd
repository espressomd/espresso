
from espresso cimport *
## Here we create something to handle particles
cimport numpy as np
from utils cimport *

cdef extern from "../src/particle_data.h":
  int place_particle(int part, double p[3])
  ctypedef struct IntList:
    pass
  ctypedef struct ParticleProperties:
    pass
  ctypedef struct ParticlePosition:
    double p[3]
  ctypedef struct ParticleLocal:
    pass
  ctypedef struct ParticleMomentum:
    pass
  ctypedef struct ParticleForce:
    pass
  ctypedef struct ParticleLocal:
    pass
  ctypedef struct Particle:
#    ParticleProperties p
    ParticlePosition r
#    ParticleMomentum m
#    ParticleForce f
#    ParticleLocal l

#    if cython.cmacro(cython.defined(LB)):
#      ParticleLatticeCoupling lc
    IntList bl
#    if cython.cmacro(cython.defined(EXCLUSIONS)):
#    IntList el
  int get_particle_data(int part, Particle *data)


cdef class ParticleHandle:
  cdef public int id
  cdef bint valid
  cdef Particle particleData
  cdef int update_particle_data(self)
  



