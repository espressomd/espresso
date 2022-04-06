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
from .utils cimport Vector4d, Vector3d, Vector3i, Span, Quaternion
from libcpp cimport bool
from libcpp.vector cimport vector  # import std::vector as vector

include "myconfig.pxi"

from .utils cimport Span


cdef extern from "Particle.hpp":
    cppclass BondView:
        int bond_id()
        Span[const int] partner_ids()

    ctypedef struct particle_parameters_swimming "ParticleParametersSwimming":
        bool swimming
        double f_swim
        double v_swim
        int push_pull
        double dipole_length

    # Note: Conditional compilation is not possible within ctypedef blocks.
    # Therefore, only member variables are imported here, which are always compiled into ESPResSo.
    # For all other properties, getter-functions have to be used on the c
    # level.        

    ctypedef struct particle "Particle":
        Vector3d calc_dip()
        int type()
        int identity()
        double mass()
        int mol_id()
        Vector3d pos()
        Vector3d calc_director()
        Vector3d force()
        Vector3d v()
        Vector3i image_box()
        double lees_edwards_offset()
        int lees_edwards_flag()
        Vector3d rinertia()
        Vector3d mu_E()
        double q()
        Vector3d omega()
        Vector3d torque()
        Quaternion[double] quat()
        double dipm()
        bint is_virtual()
        Vector3d ext_force()
        Vector3d ext_torque()
        vector[int] exclusions_as_vector() except +
        bool has_exclusion(int pid) except +
        particle_parameters_swimming swimming()

cdef extern from "particle_data.hpp":
    # Setter/getter/modifier functions functions
    void prefetch_particle_data(vector[int] ids)

    void set_particle_v(int part, const Vector3d & v)

    void set_particle_f(int part, const Vector3d & f)
    void set_particle_lees_edwards_offset(int, const double)

    IF MASS:
        void set_particle_mass(int part, double mass)

    IF ROTATIONAL_INERTIA:
        void set_particle_rotational_inertia(int part, const Vector3d & rinertia)

    IF ROTATION:
        void set_particle_rotation(int part, const Vector3i & flag)
        Vector3i get_particle_rotation(const particle * p)

    void set_particle_q(int part, double q)

    IF LB_ELECTROHYDRODYNAMICS:
        void set_particle_mu_E(int part, const Vector3d & mu_E)

    void set_particle_type(int part, int type)

    void set_particle_mol_id(int part, int mid)

    IF ROTATION:
        void set_particle_quat(int part, const Quaternion[double] & quat)
        void set_particle_director(int part, const Vector3d & director)
        void set_particle_omega_lab(int part, const Vector3d & omega)
        void set_particle_omega_body(int part, const Vector3d & omega)
        void set_particle_torque_lab(int part, const Vector3d & torque)

    IF DIPOLES:
        void set_particle_dip(int part, const Vector3d & dip)
        void set_particle_dipm(int part, double dipm)

    IF VIRTUAL_SITES:
        void set_particle_virtual(int part, int isVirtual)

    IF THERMOSTAT_PER_PARTICLE:
        IF PARTICLE_ANISOTROPY:
            void set_particle_gamma(int part, const Vector3d & gamma)
            Vector3d get_particle_gamma(const particle * p)
        ELSE:
            void set_particle_gamma(int part, double gamma)
            double get_particle_gamma(const particle * p)

        IF ROTATION:
            IF PARTICLE_ANISOTROPY:
                void set_particle_gamma_rot(int part, const Vector3d & gamma_rot)
                Vector3d get_particle_gamma_rot(const particle * p)
            ELSE:
                void set_particle_gamma_rot(int part, double gamma)
                double get_particle_gamma_rot(const particle * p)

    IF VIRTUAL_SITES_RELATIVE:
        Quaternion[double] get_particle_vs_relative(const particle * p, int & vs_relative_to, double & vs_distance)
        Quaternion[double] get_particle_vs_quat(const particle * p)
        void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance, const Quaternion[double] & rel_ori)
        void set_particle_vs_quat(int part, const Quaternion[double] & vs_quat)

    IF EXTERNAL_FORCES:
        IF ROTATION:
            void set_particle_ext_torque(int part, const Vector3d & torque)

        void set_particle_ext_force(int part, const Vector3d & force)

        void set_particle_fix(int part, const Vector3i & flag)
        Vector3i get_particle_fix(const particle * p)

    void delete_particle_bond(int part, Span[const int] bond)
    void delete_particle_bonds(int part)
    void add_particle_bond(int part, Span[const int] bond)
    const vector[BondView] & get_particle_bonds(int part)

    IF EXCLUSIONS:
        int change_exclusion(int part, int part2, int _delete)

    IF ENGINE:
        void set_particle_swimming(int part, particle_parameters_swimming swim)

    void remove_all_bonds_to(int part)

cdef extern from "particle_node.hpp":
    void place_particle(int p_id, const Vector3d & pos) except +

    void remove_particle(int part) except +

    void remove_all_particles() except +

    bool particle_exists(int part)

    int get_particle_node(int pid) except +

    const particle & get_particle_data(int pid) except +

    vector[int] get_particle_ids() except +

    int get_maximal_particle_id()
    int get_n_part()

cdef extern from "virtual_sites.hpp":
    IF VIRTUAL_SITES_RELATIVE == 1:
        void vs_relate_to(int part_num, int relate_to) except +

cdef extern from "rotation.hpp":
    Vector3d convert_vector_body_to_space(const particle & p, const Vector3d & v)
    Vector3d convert_vector_space_to_body(const particle & p, const Vector3d & v)
    void rotate_particle(int pid, const Vector3d & axis, double angle)

cdef extern from "bonded_interactions/rigid_bond.hpp":
    extern int n_rigidbonds

cdef class ParticleHandle:
    cdef public int _id
    cdef const particle * particle_data
    cdef int update_particle_data(self) except -1

cdef class _ParticleSliceImpl:
    cdef public id_selection
    cdef int _chunk_size
