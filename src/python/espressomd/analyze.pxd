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

# For C-extern Analysis

from .utils cimport Vector3i, Vector3d, Vector9d, Span
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp cimport bool as cbool

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass array4 "std::array<double, 4>":
        array4() except+
        double & operator[](size_t)

    cdef cppclass array2 "std::array<double, 2>":
        array2() except+
        double & operator[](size_t)

cdef extern from "PartCfg.hpp":
    cppclass PartCfg:
        pass

cdef extern from "partCfg_global.hpp":
    PartCfg & partCfg()

cdef extern from "particle_data.hpp":
    int max_seen_particle_type

cdef extern from "Observable_stat.hpp":
    cdef cppclass Observable_stat_wrapper:
        cbool v_comp
        Span[double] kinetic
        Span[double] bonded
        Span[double] non_bonded
        Span[double] coulomb
        Span[double] dipolar
        Span[double] virtual_sites
        Span[double] external_fields
        double accumulate(...)
        double accumulate2 "accumulate"(double acc, size_t column)
        double accumulate1 "accumulate"(double acc)
        Span[double] bonded_contribution(int bond_id)
        Span[double] non_bonded_contribution(int type1, int type2)

    cdef cppclass Observable_stat_non_bonded_wrapper:
        Span[double] non_bonded_intra_contribution(int type1, int type2)
        Span[double] non_bonded_inter_contribution(int type1, int type2)

cdef extern from "statistics.hpp":
    cdef vector[double] calc_structurefactor(PartCfg & , const vector[int] & p_types, int order)
    cdef vector[vector[double]] modify_stucturefactor(int order, double * sf)
    cdef double mindist(PartCfg & , const vector[int] & set1, const vector[int] & set2)
    cdef vector[int] nbhood(PartCfg & , const Vector3d & pos, double r_catch, const Vector3i & planedims)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef vector[double] centerofmass(PartCfg & , int part_type)

    Vector3d angularmomentum(PartCfg & , int p_type)

    void momentofinertiamatrix(PartCfg & , int p_type, double * MofImatrix)

    void calc_part_distribution(
        PartCfg &, const vector[int] & p1_types, const vector[int] & p2_types,
        double r_min, double r_max, int r_bins, bint log_flag, double * low,
        double * dist)

cdef extern from "statistics_chain.hpp":
    array4 calc_re(int, int, int)
    array4 calc_rg(int, int, int) except +
    array2 calc_rh(int, int, int)

cdef extern from "pressure_inline.hpp":
    cdef Observable_stat_wrapper obs_scalar_pressure
    cdef Observable_stat_non_bonded_wrapper obs_scalar_pressure_non_bonded
    cdef Observable_stat_wrapper obs_pressure_tensor
    cdef Observable_stat_non_bonded_wrapper obs_pressure_tensor_non_bonded

cdef extern from "pressure.hpp":
    cdef void update_pressure(cbool)

cdef extern from "energy_inline.hpp":
    cdef Observable_stat_wrapper obs_energy

cdef extern from "energy.hpp":
    cdef void update_energy()
    double calculate_current_potential_energy_of_system()

cdef extern from "dpd.hpp":
    Vector9d dpd_stress()
