from espresso cimport *
## Here we create something to handle particles
cimport numpy as np
from utils cimport *
include "myconfig.pxi"

cdef extern from "particle_data.hpp":
  int place_particle(int part, double p[3])
  ctypedef struct IntList:
    pass
  ctypedef struct ParticleProperties:
    int type
    double q
    pass
  ctypedef struct ParticlePosition:
    double p[3]
  ctypedef struct ParticleLocal:
    pass
  ctypedef struct ParticleMomentum:
    double v[3]
    pass
  ctypedef struct ParticleForce:
    double f[3]
    pass
  ctypedef struct ParticleLocal:
    pass
  ctypedef struct Particle:
    ParticleProperties p
    ParticlePosition r
    ParticleMomentum m
    ParticleForce f
    ParticleLocal l

#    if cython.cmacro(cython.defined(LB)):
#      ParticleLatticeCoupling lc
    IntList bl
#    if cython.cmacro(cython.defined(EXCLUSIONS)):
#    IntList el
  int get_particle_data(int part, Particle *data)
  int set_particle_type(int part, int type)
  int set_particle_v(int part, double v[3])
  int set_particle_f(int part, double F[3])
  IF ELECTROSTATICS == 1:
    int set_particle_q(int part, double q)
    

cdef class ParticleHandle:
  cdef public int id
  cdef bint valid
  cdef Particle particleData
  cdef int update_particle_data(self)

