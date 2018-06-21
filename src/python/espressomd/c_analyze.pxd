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
# For C-extern Analysis

#
# TODO: Merge blocks of cdef extern for same headers.
#

from __future__ import print_function, absolute_import
cimport numpy as np
from espressomd.utils cimport *
from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp.map cimport map  # import std::map as map

cdef extern from "PartCfg.hpp":
    cppclass PartCfg:
        pass

cdef extern from "partCfg_global.hpp":
    PartCfg & partCfg()

cdef extern from "particle_data.hpp":
    int max_seen_particle_type

cdef extern from "statistics.hpp":
    cdef void calc_structurefactor(PartCfg &, int * p_types, int n_types, int order, double ** sf)
    cdef vector[vector[double]] modify_stucturefactor(int order, double * sf)

cdef extern from "statistics.hpp":
    ctypedef struct Observable_stat:
        int init_status
        double_list data
        int n_coulomb
        int n_dipolar
        int n_non_bonded
        int n_virtual_sites
        double * bonded
        double * non_bonded
        double * coulomb
        double * dipolar
        double * virtual_sites

cdef extern from "statistics.hpp":
    ctypedef struct Observable_stat_non_bonded:
        pass
    cdef double mindist(PartCfg &, const int_list & set1, const int_list & set2)
    cdef double min_distance2(double pos1[3], double pos2[3])
    cdef int_list nbhood(PartCfg &, double pos[3], double r_catch, int planedims[3])
    cdef double distto(PartCfg &, double pos[3], int pid)
    cdef double * obsstat_bonded(Observable_stat * stat, int j)
    cdef double * obsstat_nonbonded(Observable_stat * stat, int i, int j)
    cdef double * obsstat_nonbonded_inter(Observable_stat_non_bonded * stat, int i, int j)
    cdef double * obsstat_nonbonded_intra(Observable_stat_non_bonded * stat, int i, int j)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef vector[double] centerofmass(PartCfg &, int part_type)
    cdef int calc_cylindrical_average(PartCfg &, vector[double] center, vector[double] direction, double length,
                                      double radius, int bins_axial, int bins_radial, vector[int] types,
                                      map[string, vector[vector[vector[double]]]] & distribution)

cdef extern from "pressure.hpp":
    cdef Observable_stat total_pressure
    cdef Observable_stat_non_bonded total_pressure_non_bonded
    cdef Observable_stat total_p_tensor
    cdef Observable_stat_non_bonded total_p_tensor_non_bonded
    cdef void update_pressure(int)
    cdef int analyze_local_stress_tensor(int * periodic, double * range_start, double * range, int * bins, double_list * local_stress_tensor)

cdef extern from "energy.hpp":
    cdef Observable_stat total_energy
    cdef Observable_stat_non_bonded total_energy_non_bonded
    cdef void master_energy_calc()
    cdef void init_energies(Observable_stat * stat)

cdef extern from "statistics_chain.hpp":
    int chain_start
    int chain_n_chains
    int chain_length
    void calc_re(PartCfg&, double ** re)
    void calc_rg(PartCfg&, double ** rg)
    void calc_rh(PartCfg&, double ** rh)

cdef extern from "interaction_data.hpp":
    int n_bonded_ia

cdef extern from "statistics.hpp":
    void calc_rdf(PartCfg &, vector[int] p1_types, vector[int] p2_types,
                  double r_min, double r_max, int r_bins, vector[double] rdf)

    void calc_rdf_av(PartCfg &, vector[int] p1_types, vector[int] p2_types,
                     double r_min, double r_max, int r_bins, vector[double] rdf, int n_conf)

    void angularmomentum(PartCfg &, int p_type, double * com)
    void momentofinertiamatrix(PartCfg &, int p_type, double * MofImatrix)
    void analyze_rdfchain(PartCfg &, double r_min, double r_max, int r_bins, double ** f1, double ** f2, double ** f3)

cdef extern from "statistics.hpp":
    int n_part
    int n_part_conf
    int n_configs
    void analyze_append(PartCfg &)

cdef extern from "statistics.hpp":
    void calc_part_distribution(PartCfg &, int *p1_types, int n_p1, int *p2_types, int n_p2,
                                double r_min, double r_max, int r_bins, int log_flag, 
                                double *low, double *dist)
