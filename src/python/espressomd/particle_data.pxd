#
# Copyright (C) 2013-2019 The ESPResSo project
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
# Here we create something to handle particles
from .utils cimport Vector4d, Vector3d, Vector3i, Span
from libcpp cimport bool
from libcpp.vector cimport vector  # import std::vector as vector
from libc cimport stdint

include "myconfig.pxi"

from .utils cimport Span

# Import particle data structures and setter functions from particle_data.hpp
cdef extern from "particle_data.hpp":
    cppclass BondView:
        int bond_id()
        Span[const int] partner_ids()

    # DATA STRUCTURES
    stdint.uint8_t ROTATION_X
    stdint.uint8_t ROTATION_Y
    stdint.uint8_t ROTATION_Z

    # Note: Conditional compilation is not possible within ctypedef blocks.
    # Therefore, only member variables are imported here, which are always compiled into ESPResSo.
    # For all other properties, getter-functions have to be used on the c
    # level.
    ctypedef struct particle_properties "ParticleProperties":
        int    identity
        int    mol_id
        int    type
        double mass
        stdint.uint8_t rotation

    ctypedef struct particle_position "ParticlePosition":
        Vector3d p
        Vector3d calc_director()

    ctypedef struct particle_force "ParticleForce":
        Vector3d f

    ctypedef struct particle_momentum "ParticleMomentum":
        Vector3d v

    ctypedef struct particle_local "ParticleLocal":
        Vector3i i

    ctypedef struct particle "Particle":
        particle_properties p
        particle_position r
        particle_momentum m
        particle_force f
        particle_local l
        vector[int] exclusions() except +
        Vector3d calc_dip()

    IF ENGINE:
        ctypedef struct particle_parameters_swimming "ParticleParametersSwimming":
            bool swimming
            double f_swim
            double v_swim
            int push_pull
            double dipole_length

    # Setter/getter/modifier functions functions
    void prefetch_particle_data(vector[int] ids)

    int place_particle(int part, double p[3])

    void set_particle_v(int part, double v[3])

    void set_particle_f(int part, const Vector3d & F)

    IF ROTATION:
        void set_particle_rotation(int part, int rot)

    IF MASS:
        void set_particle_mass(int part, double mass)

    IF ROTATIONAL_INERTIA:
        void set_particle_rotational_inertia(int part, double rinertia[3])
        void pointer_to_rotational_inertia(const particle * p, const double * & res)

    IF ROTATION:
        void set_particle_rotation(int part, int rot)

    void set_particle_q(int part, double q)

    IF LB_ELECTROHYDRODYNAMICS:
        void set_particle_mu_E(int part, const Vector3d & mu_E)
        void get_particle_mu_E(int part, Vector3d & mu_E)

    void set_particle_type(int part, int type)

    void set_particle_mol_id(int part, int mid)

    IF ROTATION:
        void set_particle_quat(int part, double quat[4])
        void pointer_to_quat(const particle * p, const double * & res)
        void set_particle_omega_lab(int part, Vector3d omega)
        void set_particle_omega_body(int part, Vector3d omega)
        void set_particle_torque_lab(int part, Vector3d torque)
        void pointer_to_omega_body(const particle * p, const double * & res)
        Vector3d get_torque_body(const particle p)

    IF DIPOLES:
        void set_particle_dip(int part, double dip[3])

        void set_particle_dipm(int part, double dipm)
        void pointer_to_dipm(const particle * P, const double * & res)

    IF VIRTUAL_SITES:
        void set_particle_virtual(int part, int isVirtual)
        void pointer_to_virtual(const particle * P, const bint * & res)

    IF LANGEVIN_PER_PARTICLE:
        void set_particle_temperature(int part, double T)
        void pointer_to_temperature(const particle * p, const double * & res)

        IF PARTICLE_ANISOTROPY:
            void set_particle_gamma(int part, Vector3d gamma)
        ELSE:
            void set_particle_gamma(int part, double gamma)

        void pointer_to_gamma(const particle * p, const double * & res)

        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                void set_particle_gamma_rot(int part, Vector3d gamma_rot)
            ELSE:
                void set_particle_gamma_rot(int part, double gamma)

            void pointer_to_gamma_rot(const particle * p, const double * & res)

    IF VIRTUAL_SITES_RELATIVE:
        void pointer_to_vs_relative(const particle * P, const int * & res1, const double * & res2, const double * & res3)
        void pointer_to_vs_quat(const particle * P, const double * & res)
        void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance, Vector4d rel_ori)
        void set_particle_vs_quat(int part, Vector4d vs_quat)

    void pointer_to_q(const particle * P, const double * & res)

    IF EXTERNAL_FORCES:
        IF ROTATION:
            void set_particle_ext_torque(int part, const Vector3d & torque)
            void pointer_to_ext_torque(const particle * P, const double * & res2)

        void set_particle_ext_force(int part, const Vector3d & force)
        void pointer_to_ext_force(const particle * P, const double * & res2)

        void set_particle_fix(int part, stdint.uint8_t flag)
        void pointer_to_fix(const particle * P, const stdint.uint8_t * & res)

    void delete_particle_bond(int part, Span[const int] bond)
    void delete_particle_bonds(int part)
    void add_particle_bond(int part, Span[const int] bond)
    const vector[BondView] & get_particle_bonds(int part)

    IF EXCLUSIONS:
        int change_exclusion(int part, int part2, int _delete)
        void remove_all_exclusions()

    IF ENGINE:
        void set_particle_swimming(int part, particle_parameters_swimming swim)
        void pointer_to_swimming(const particle * p, const particle_parameters_swimming * & swim)

    int remove_particle(int part) except +

    void remove_all_particles() except +

    void remove_all_bonds_to(int part)

    bool particle_exists(int part)

    int get_particle_node(int id) except +

    const particle & get_particle_data(int id) except +

    vector[int] get_particle_ids() except +

    int get_maximal_particle_id()
    int get_n_part()

cdef extern from "virtual_sites.hpp":
    IF VIRTUAL_SITES_RELATIVE == 1:
        void vs_relate_to(int part_num, int relate_to)

cdef extern from "rotation.hpp":
    Vector3d convert_vector_body_to_space(const particle & p, const Vector3d & v)
    Vector3d convert_vector_space_to_body(const particle & p, const Vector3d & v)
    void rotate_particle(int id, Vector3d axis, double angle)

cdef class ParticleHandle:
    cdef public int _id
    cdef const particle * particle_data
    cdef int update_particle_data(self) except -1

cdef class _ParticleSliceImpl:
    cdef public id_selection
    cdef int _chunk_size
