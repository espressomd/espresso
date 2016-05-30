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
from _system cimport *
# Here we create something to handle particles
cimport numpy as np
from utils cimport *
from libcpp cimport bool

include "myconfig.pxi"

# Import particle data structures and setter functions from particle_data.hpp

cdef extern from "particle_data.hpp":

    # DATA STRUCTURES

    # Note: Conditional compilation is not possible within ctypedef blocks.
    # Therefore, only member variables are imported here, which are always compiled into Espresso.
    # For all other properties, getter-funcionts have to be used on the c
    # level.

    ctypedef struct particle_properties "ParticleProperties":
        int    identity
        int    mol_id
        int    type

    ctypedef struct particle_position "ParticlePosition":
        double p[3]

    ctypedef struct particle_force "ParticleForce":
        double f[3]

    ctypedef struct particle_momentum "ParticleMomentum":
        double v[3]

    ctypedef struct particle_local "ParticleLocal":
        pass

    ctypedef struct particle "Particle":
        particle_properties p
        particle_position r
        particle_momentum m
        particle_force f
        particle_local l
        int_list bl

    IF ENGINE:
        IF LB or LB_GPU:
            ctypedef struct particle_parameters_swimming "ParticleParametersSwimming":
                bool swimming
                double f_swim
                double v_swim
                int push_pull
                double dipole_length
                double v_center[3]
                double v_source[3]
                double rotational_friction
        ELSE:
            ctypedef struct particle_parameters_swimming "ParticleParametersSwimming":
                bool swimming
                double f_swim
                double v_swim

    # Setter/getter/modifier functions functions

    int get_particle_data(int part, particle * data)

    int place_particle(int part, double p[3])

    int set_particle_v(int part, double v[3])

    int set_particle_f(int part, double F[3])

    int set_particle_mass(int part, double mass)

    int set_particle_solvation(int part, double * solvation)

    IF ROTATION_PER_PARTICLE == 1:
        int set_particle_rotation(int part, int rot)

    IF MULTI_TIMESTEP:
        int set_particle_smaller_timestep(int part, int small_timestep)
        void pointer_to_smaller_timestep(particle * p, int * & res)

    IF MASS:
        int set_particle_mass(int part, double mass)
        void pointer_to_mass(particle * p, double * & res)

    IF SHANCHEN:
        int set_particle_solvation(int part, double * solvation)
        void pointer_to_solvation(particle * p, double * & res)

    IF ROTATIONAL_INERTIA:
        int set_particle_rotational_inertia(int part, double rinertia[3])
        void pointer_to_rotational_inertia(particle * p, double * & res)

    IF ROTATION_PER_PARTICLE:
        int set_particle_rotation(int part, int rot)
        void pointer_to_rotation(particle * p, short int * & res)

    IF ELECTROSTATICS:
        int set_particle_q(int part, double q)

    int set_particle_mu_E(int part, double mu_E[3])

    int set_particle_type(int part, int type)

    int set_particle_mol_id(int part, int mid)

    IF ROTATION:
        int set_particle_quat(int part, double quat[4])
        void pointer_to_quat(particle * p, double * & res)
        void pointer_to_quatu(particle * p, double * & res)
        int set_particle_omega_lab(int part, double omega[3])
        int set_particle_omega_body(int part, double omega[3])
        int set_particle_torque_lab(int part, double torque[3])
        int set_particle_torque_body(int part, double torque[3])
        void pointer_to_omega_body(particle * p, double * & res)
        void pointer_to_torque_lab(particle * p, double * & res)

    IF MASS == 1:
        void pointer_to_mass(particle * p, double * & res)

    IF DIPOLES:
        int set_particle_dip(int part, double dip[3])
        void pointer_to_dip(particle * P, double * & res)

        int set_particle_dipm(int part, double dipm)
        void pointer_to_dipm(particle * P, double * & res)

    IF VIRTUAL_SITES:
        int set_particle_virtual(int part, int isVirtual)
        void pointer_to_virtual(particle * P, int * & res)

    IF LANGEVIN_PER_PARTICLE:
        int set_particle_temperature(int part, double T)
        void pointer_to_temperature(particle * p, double * & res)

        int set_particle_gamma(int part, double gamma)
        void pointer_to_gamma(particle * p, double * & res)

    IF VIRTUAL_SITES_RELATIVE:
        void pointer_to_vs_relative(particle * P, int * & res1, double * & res2, double * & res3)

    IF ELECTROSTATICS:
        void pointer_to_q(particle * P, double * & res)

    IF EXTERNAL_FORCES:
        IF ROTATION:
            int set_particle_ext_torque(int part, int flag, double torque[3])
            void pointer_to_ext_torque(particle * P, int * & res1, double * & res2)

        int set_particle_ext_force(int part, int flag, double force[3])
        void pointer_to_ext_force(particle * P, int * & res1, double * & res2)

        int set_particle_fix(int part,  int flag)
        void pointer_to_fix(particle * P, int * & res)

    int change_particle_bond(int part, int * bond, int _delete)

    IF EXCLUSIONS:
        int change_exclusion(int part, int part2, int _delete)
        void pointer_to_exclusions(particle * p, int * & res1, int * & res2)

        void remove_all_exclusions()

    IF ENGINE:
        int set_particle_swimming(int part, particle_parameters_swimming swim)
        void pointer_to_swimming(particle * p, particle_parameters_swimming * & swim)

    int remove_particle(int part)

    void remove_all_particles()

    void remove_all_bonds_to(int part)

    bool particle_exists(int part)

cdef extern from "virtual_sites_relative.hpp":
    IF VIRTUAL_SITES_RELATIVE == 1:
        int vs_relate_to(int part_num, int relate_to)
        int set_particle_vs_relative(int part, int vs_relative_to, double vs_distance, double * vs_quat)

cdef extern from "rotation.hpp":
    void convert_omega_body_to_space(particle * p, double * omega)
    void convert_torques_body_to_space(particle * p, double * torque)

# The bonded_ia_params stuff has to be included here, because the setter/getter
# of the particles' bond property needs to now about the correct number of
# bond partners
cdef extern from "interaction_data.hpp":
    ctypedef struct bonded_ia_parameters "Bonded_ia_parameters":
        int num
        pass
    bonded_ia_parameters * bonded_ia_params
    cdef int n_bonded_ia

cdef class ParticleHandle(object):

    cdef public int id
    cdef bint valid
    cdef particle particle_data
    cdef int update_particle_data(self) except -1


cdef class ParticleSlice:

    cdef particle particle_data
    cdef int update_particle_data(self,id) except -1
    cdef public id_selection
