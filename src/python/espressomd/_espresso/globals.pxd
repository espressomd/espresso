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
include "myconfig.pxi"

cdef extern from "global.hpp":
    int FIELD_MAXNUMCELLS
    int FIELD_MINNUMCELLS
    int FIELD_NODEGRID
    int FIELD_NPTISO_PISTON
    int FIELD_NPTISO_PDIFF
    int FIELD_PERIODIC
    int FIELD_SIMTIME

cdef extern from "communication.hpp":
    extern int n_nodes
    void mpi_set_smaller_time_step(double smaller_time_step)
    void mpi_set_time_step(double time_step)
    void mpi_bcast_parameter(int p)

cdef extern from "integrate.hpp":
    double time_step
    extern int integ_switch
    extern double sim_time
    extern double smaller_time_step
    extern double verlet_reuse

cdef extern from "verlet.hpp":
    double skin

cdef extern from "lattice.hpp":
    extern int lattice_switch

cdef extern from "domain_decomposition.hpp":
    ctypedef struct IA_Neighbor:
        pass
    ctypedef struct IA_Neighbor_List:
        pass
    ctypedef struct  DomainDecomposition:
        int use_vList
        int cell_grid[3]
        double cell_size[3]

    extern DomainDecomposition dd
    extern int max_num_cells
    extern int min_num_cells
    extern double max_skin
    int calc_processor_min_num_cells()


cdef extern from "particle_data.hpp":
    extern int n_part


cdef extern from "interaction_data.hpp":
    double dpd_gamma
    double dpd_r_cut
    extern double max_cut
    extern int max_seen_particle
    extern int n_particle_types
    extern double max_cut_nonbonded
    extern double max_cut_bonded


cdef extern from "thermostat.hpp":
    double langevin_gamma
    extern double nptiso_gamma0
    extern double nptiso_gammav
    extern double temperature
    extern int thermo_switch


cdef extern from "dpd.hpp":
    extern int dpd_wf
    extern double dpd_tgamma
    extern double dpd_tr_cut
    extern int dpd_twf


IF LB:
    cdef extern from "lb.hpp":
        ctypedef struct LB_Parameters:
            double tau
        extern LB_Parameters lbpar

IF LB_GPU:
    cdef extern from "lbgpu.hpp":
        ctypedef struct LB_parameters_gpu:
            double tau
        extern LB_parameters_gpu lbpar_gpu

cdef extern from "cells.hpp":
    extern double max_range
    ctypedef struct CellStructure:
        int type
    CellStructure cell_structure

cdef extern from "layered.hpp":
    extern int n_layers

cdef extern from "rattle.hpp":
    extern int n_rigidbonds


cdef extern from "tuning.hpp":
    extern int timing_samples

cdef extern from "imd.hpp":
    extern int transfer_rate


cdef extern from "grid.hpp":
    double box_l[3]
    double local_box_l[3]
    extern int node_grid[3]
    extern int periodic
    extern double min_box_l

cdef extern from "npt.hpp":
    ctypedef struct nptiso_struct:
        double p_ext
        double p_inst
        double p_inst_av
        double p_diff
        double piston
    extern nptiso_struct nptiso

cdef extern from "statistics.hpp":
    extern int n_configs

cdef extern from "reaction.hpp":
    ctypedef struct  reaction_struct:
        int reactant_type
        int product_type
        int catalyzer_type
        double range
        double ct_rate
        double eq_rate
        int sing_mult
        int swap

    cdef extern reaction_struct reaction
