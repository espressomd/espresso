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

from .utils cimport Vector3i, Vector3d, Vector9d
from libcpp.vector cimport vector  # import std::vector as vector

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

cdef extern from "statistics.hpp":
    int get_n_part_conf()
    int get_n_configs()

    ctypedef struct Observable_stat:
        int init_status
        vector[double] data
        int n_coulomb
        int n_dipolar
        int n_non_bonded
        int n_virtual_sites
        double * bonded
        double * non_bonded
        double * coulomb
        double * dipolar
        double * virtual_sites
        double * external_fields

    ctypedef struct Observable_stat_non_bonded:
        pass

    cdef vector[double] calc_structurefactor(PartCfg & , const vector[int] & p_types, int order)
    cdef vector[vector[double]] modify_stucturefactor(int order, double * sf)
    cdef double mindist(PartCfg & , const vector[int] & set1, const vector[int] & set2)
    cdef vector[int] nbhood(PartCfg & , const Vector3d & pos, double r_catch, const Vector3i & planedims)
    cdef double * obsstat_bonded(Observable_stat * stat, int j)
    cdef double * obsstat_nonbonded(Observable_stat * stat, int i, int j)
    cdef double * obsstat_nonbonded_inter(Observable_stat_non_bonded * stat, int i, int j)
    cdef double * obsstat_nonbonded_intra(Observable_stat_non_bonded * stat, int i, int j)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef vector[double] centerofmass(PartCfg & , int part_type)

    void calc_rdf(PartCfg & , vector[int] p1_types, vector[int] p2_types,
                  double r_min, double r_max, int r_bins, vector[double] rdf)

    void calc_rdf_av(PartCfg & , vector[int] p1_types, vector[int] p2_types,
                     double r_min, double r_max, int r_bins, vector[double] rdf,
                     int n_conf)

    Vector3d angularmomentum(PartCfg & , int p_type)

    void momentofinertiamatrix(PartCfg & , int p_type, double * MofImatrix)

    void analyze_append(PartCfg & )

    void calc_part_distribution(
        PartCfg &, const vector[int] & p1_types, const vector[int] & p2_types,
        double r_min, double r_max, int r_bins, bint log_flag, double * low,
        double * dist)

cdef extern from "statistics_chain.hpp":
    array4 calc_re(int, int, int)
    array4 calc_rg(int, int, int) except +
    array2 calc_rh(int, int, int)

cdef extern from "pressure.hpp":
    cdef Observable_stat total_pressure
    cdef Observable_stat_non_bonded total_pressure_non_bonded
    cdef Observable_stat total_p_tensor
    cdef Observable_stat_non_bonded total_p_tensor_non_bonded
    cdef void update_pressure(int)

cdef extern from "energy.hpp":
    cdef Observable_stat total_energy
    cdef Observable_stat_non_bonded total_energy_non_bonded
    cdef void master_energy_calc()
    cdef void init_energies(Observable_stat * stat)
    double calculate_current_potential_energy_of_system()

cdef extern from "dpd.hpp":
    Vector9d dpd_stress()
