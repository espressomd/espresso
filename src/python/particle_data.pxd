from espresso cimport *
## Here we create something to handle particles
cimport numpy as np
from utils cimport *
include "myconfig.pxi"


# Import particle data structures and setter functions from particle_data.hpp

cdef extern from "particle_data.hpp":
  
  # DATA STRUCTURES

  # Note: Conditional compilation is not possible within ctypedef blocks. 
  # Therefore, only member variables are imported here, which are always compiled into Espresso.
  # For all other properties, getter-funcionts have to be used on the c level.
  
  ctypedef struct ParticleProperties:
    int    identity
    int    mol_id
    int    type
    
  ctypedef struct ParticlePosition:
    double p[3]


  ctypedef struct ParticleForce:
    double f[3]

  ctypedef struct ParticleMomentum:
    double v[3]

  ctypedef struct ParticleLocal:
    pass


  ctypedef struct Particle: 
    ParticleProperties p
    ParticlePosition r
    ParticleMomentum m
    ParticleForce f
    ParticleLocal l
    IntList bl


  # Setter functions + a few other oens

  int get_particle_data(int part, Particle *data)
  
  int place_particle(int part, double p[3])
  
  int set_particle_v(int part, double v[3])
  
  
  int set_particle_f(int part, double F[3])
  
  int set_particle_mass(int part, double mass)
  
  int set_particle_solvation(int part, double* solvation)
  
  
  IF ROTATIONAL_INERTIA == 1:
    int set_particle_rotational_inertia(int part, double rinertia[3])
  
  IF ROTATION_PER_PARTICLE == 1:
    int set_particle_rotation(int part, int rot)
  
  
  int set_particle_q(int part, double q)
  
  int set_particle_mu_E(int part, double mu_E[3])
  
  int set_particle_type(int part, int type)
  
  int set_particle_mol_id(int part, int mid)
  
  IF ROTATION == 1:
    int set_particle_quat(int part, double quat[4])
  
    int set_particle_omega_lab(int part, double omega[3])
  
    int set_particle_omega_body(int part, double omega[3])
  
    int set_particle_torque_lab(int part, double torque[3])
  
    int set_particle_torque_body(int part, double torque[3])
  
  IF DIPOLES == 1:
    int set_particle_dip(int part, double dip[3])
  
    int set_particle_dipm(int part, double dipm)
  
  IF VIRTUAL_SITES == 1:
    int set_particle_virtual(int part,int isVirtual)
  
  IF LANGEVIN_PER_PARTICLE == 1:
    int set_particle_temperature(int part, double T)
  
    int set_particle_gamma(int part, double gamma)
  
  IF EXTERNAL_FORCES == 1:
    IF ROTATION == 1:
      int set_particle_ext_torque(int part, int flag, double torque[3])
    int set_particle_ext_force(int part, int flag, double force[3])
    int set_particle_fix(int part,  int flag)
  
  int change_particle_bond(int part, int *bond, int _delete)
  
  IF EXCLUSIONS == 1:
    int change_exclusion(int part, int part2, int _delete)
  
    void remove_all_exclusions()
  
  
  int remove_particle(int part)
  
  void remove_all_particles()
  
  void remove_all_bonds_to(int part)
  
  
cdef class ParticleHandle(object):
  cdef public int id
  cdef bint valid
  cdef Particle particleData
  cdef int update_particle_data(self)

