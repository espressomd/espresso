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
from __future__ import print_function, absolute_import
from espressomd.system cimport *
# Here we create something to handle particles
cimport numpy as np
from espressomd.utils cimport Vector3d,int_list
from espressomd.utils import array_locked
from libcpp cimport bool
from libcpp.memory cimport unique_ptr

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
        double mass

    ctypedef struct particle_position "ParticlePosition":
        Vector3d p

    ctypedef struct particle_force "ParticleForce":
        Vector3d f

    ctypedef struct particle_momentum "ParticleMomentum":
        Vector3d v

    ctypedef struct particle_local "ParticleLocal":
        int i[3]

    ctypedef struct particle "Particle":
        particle_properties p
        particle_position r
        particle_momentum m
        particle_force f
        particle_local l
        int_list bl
        int_list exclusions() except +

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
    const particle* get_particle_data_ptr(int part) except +
    void prefetch_particle_data(vector[int] ids)

    int place_particle(int part, double p[3])

    int set_particle_v(int part, double v[3])

    int set_particle_f(int part, const Vector3d &F)

    int set_particle_solvation(int part, double * solvation)

    IF ROTATION == 1:
        int set_particle_rotation(int part, int rot)

    IF MASS:
        int set_particle_mass(int part, double mass)

    IF SHANCHEN:
        int set_particle_solvation(int part, double * solvation)
        void pointer_to_solvation(const particle * p, const double * & res)

    IF ROTATIONAL_INERTIA:
        int set_particle_rotational_inertia(int part, double rinertia[3])
        void pointer_to_rotational_inertia(const particle * p, const double * & res)

    IF ROTATION:
        int set_particle_rotation(int part, int rot)
        void pointer_to_rotation(const particle * p, const short int * & res)

    IF ELECTROSTATICS:
        int set_particle_q(int part, double q)

    IF LB_ELECTROHYDRODYNAMICS:
        int set_particle_mu_E(int part, double mu_E[3])
        void get_particle_mu_E(int part, double (&mu_E)[3])

    int set_particle_type(int part, int type)

    int set_particle_mol_id(int part, int mid)

    IF ROTATION:
        int set_particle_quat(int part, double quat[4])
        void pointer_to_quat(const particle * p, const double * & res)
        void pointer_to_quatu(const particle * p, const double * & res)
        int set_particle_omega_lab(int part, double omega[3])
        int set_particle_omega_body(int part, double omega[3])
        int set_particle_torque_lab(int part, double torque[3])
        int set_particle_torque_body(int part, double torque[3])
        void pointer_to_omega_body(const particle * p, const double * & res)
        void pointer_to_torque_lab(const particle * p, const double * & res)

    IF MEMBRANE_COLLISION:
        int set_particle_out_direction(int part, double out_direction[3])
        void pointer_to_out_direction(particle*  p, double*& res)
        
    IF AFFINITY:
        int set_particle_affinity(int part, double bond_site[3])
        void pointer_to_bond_site(particle*  p, double*& res)

    IF MASS == 1:
        void pointer_to_mass(particle * p, double * & res)

    IF DIPOLES:
        int set_particle_dip(int part, double dip[3])
        void pointer_to_dip(const particle * P, const double * & res)

        int set_particle_dipm(int part, double dipm)
        void pointer_to_dipm(const particle * P, const double * & res)

    IF VIRTUAL_SITES:
        int set_particle_virtual(int part, int isVirtual)
        void pointer_to_virtual(const particle * P, const int * & res)

    IF LANGEVIN_PER_PARTICLE:
        int set_particle_temperature(int part, double T)
        void pointer_to_temperature(const particle * p, const double * & res)

        IF PARTICLE_ANISOTROPY:
            int set_particle_gamma(int part, Vector3d gamma)
        ELSE:
            int set_particle_gamma(int part, double gamma)

        void pointer_to_gamma(const particle * p, const double * & res)

        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                int set_particle_gamma_rot(int part, Vector3d gamma_rot)
            ELSE:
                int set_particle_gamma_rot(int part, double gamma)

            void pointer_to_gamma_rot(const particle * p, const double * & res)

    IF VIRTUAL_SITES_RELATIVE:
        void pointer_to_vs_relative(const particle * P, const int * & res1, const double * & res2, const double * & res3)
        void pointer_to_vs_quat(const particle *P, const double *& res)
        int set_particle_vs_relative(int part, int vs_relative_to, double vs_distance, double *rel_ori)
        void set_particle_vs_quat(int part, double *vs_quat)

    IF ELECTROSTATICS:
        void pointer_to_q(const particle * P, const double * & res)

    IF EXTERNAL_FORCES:
        IF ROTATION:
            int set_particle_ext_torque(int part, int flag, double torque[3])
            void pointer_to_ext_torque(const particle * P, const int * & res1, const double * & res2)

        int set_particle_ext_force(int part, int flag, double force[3])
        void pointer_to_ext_force(const particle * P, const int * & res1, const double * & res2)

        int set_particle_fix(int part,  int flag)
        void pointer_to_fix(const particle * P, const int * & res)

    int change_particle_bond(int part, int * bond, int _delete)

    IF EXCLUSIONS:
        int change_exclusion(int part, int part2, int _delete)
        void remove_all_exclusions()

    IF ENGINE:
        int set_particle_swimming(int part, particle_parameters_swimming swim)
        void pointer_to_swimming(const particle * p, const particle_parameters_swimming * & swim)

    int remove_particle(int part)

    void remove_all_particles()

    void remove_all_bonds_to(int part)

    bool particle_exists(int part)

    int get_particle_node(int id)

cdef extern from "virtual_sites.hpp":
    IF VIRTUAL_SITES_RELATIVE == 1:
        int vs_relate_to(int part_num, int relate_to)

cdef extern from "rotation.hpp":
    void convert_omega_body_to_space(const particle * p, double * omega)
    void convert_torques_body_to_space(const particle * p, double * torque)
    Vector3d convert_vector_body_to_space(const particle& p,const Vector3d& v)
    Vector3d convert_vector_space_to_body(const particle& p,const Vector3d& v)
    void rotate_particle(int id, double* axis, double angle)

cdef class ParticleHandle(object):
    cdef public int _id
    cdef const particle * particle_data
    cdef int update_particle_data(self) except -1

cdef class _ParticleSliceImpl:
    cdef public id_selection
    cdef int _chunk_size

cdef extern from "grid.hpp":
    Vector3d folded_position(const particle *)
    Vector3d unfolded_position(const particle *)
    cdef void fold_position(double *, int*)
    void unfold_position(double pos[3], int image_box[3])

cdef make_array_locked(const Vector3d &v)

