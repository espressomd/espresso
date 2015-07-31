#
# Copyright (C) 2013,2014 The ESPResSo project
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

cimport numpy as np
from utils cimport *
from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp.map cimport map  # import std::map as map

cdef extern from "particle_data.hpp":
    cdef int updatePartCfg(int bonds_flag)
    int n_particle_types

cdef extern from "statistics.hpp":
    cdef void calc_structurefactor(int type, int order, double ** sf)
    cdef vector[vector[double]] modify_stucturefactor(int order, double * sf)

cdef extern from "statistics.hpp":
    ctypedef struct Observable_stat:
        int init_status
        DoubleList data
        int n_coulomb
        int n_dipolar
        int n_non_bonded
        double * bonded
        double * non_bonded
        double * coulomb
        double * dipolar
        double * vs_relative

cdef extern from "statistics.hpp":
    ctypedef struct Observable_stat_non_bonded:
        pass
    cdef double mindist(IntList * set1, IntList * set2)
    cdef void nbhood(double pos[3], double r_catch, IntList * il, int planedims[3])
    cdef double distto(double pos[3], int pid)
    cdef double * obsstat_bonded(Observable_stat * stat, int j)
    cdef double * obsstat_nonbonded(Observable_stat * stat, int i, int j)
    cdef double * obsstat_nonbonded_inter(Observable_stat_non_bonded * stat, int i, int j)
    cdef double * obsstat_nonbonded_intra(Observable_stat_non_bonded * stat, int i, int j)
    cdef double mindist(IntList * set1, IntList * set2)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef int calc_cylindrical_average(vector[double] center, vector[double] direction, double length,
                                      double radius, int bins_axial, int bins_radial, vector[int] types,
                                      map[string, vector[vector[vector[double]]]] & distribution)


cdef extern from "pressure.hpp":
    cdef Observable_stat total_pressure
    cdef Observable_stat_non_bonded total_pressure_non_bonded
    cdef Observable_stat total_p_tensor
    cdef Observable_stat_non_bonded total_p_tensor_non_bonded
    cdef void update_pressure(int)
    cdef void analyze_pressure_all(vector[string] & pressure_labels, vector[double] & pressures, int v_comp)
    cdef double analyze_pressure(string pressure_to_calc, int v_comp)
    cdef double analyze_pressure_pair(string pressure_to_calc, int type1, int type2, int v_comp)
    cdef double analyze_pressure_single(string pressure_to_calc, int bond_or_type, int v_comp)
    cdef void analyze_stress_tensor_all(vector[string] & stressTensorLabel, vector[double] & stressTensorValues, int v_comp)
    cdef int analyze_stress_tensor(string pressure_to_calc, int v_comp, vector[double] & stress)
    cdef int analyze_stress_pair(string pressure_to_calc, int type1, int type2, int v_comp, vector[double] & stress)
    cdef int analyze_stress_single(string pressure_to_calc, int bond_or_type, int v_comp, vector[double] & stress)
    cdef int analyze_local_stress_tensor(int * periodic, double * range_start, double * range, int * bins, DoubleList * local_stress_tensor)

cdef extern from "energy.hpp":
    cdef Observable_stat total_energy
    cdef Observable_stat_non_bonded total_energy_non_bonded

cdef extern from "energy.hpp":
    cdef void master_energy_calc()
    cdef void init_energies(Observable_stat * stat)

cdef extern from "statistics_chain.hpp":
    int sortPartCfg()
    int chain_start
    int chain_n_chains
    int chain_length
    void calc_re(double ** re)
    void calc_rg(double ** rg)
    void calc_rh(double ** rh)

cdef extern from "interaction_data.hpp":
    int n_bonded_ia
