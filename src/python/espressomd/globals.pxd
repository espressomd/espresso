#
# Copyright (C) 2013-2018 The ESPResSo project
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
from libcpp cimport bool
from interactions cimport ImmersedBoundaries
from utils cimport Vector3i

cdef extern from "global.hpp":
    int FIELD_BOXL
    int FIELD_SKIN
    int FIELD_NODEGRID
    int FIELD_MAXNUMCELLS
    int FIELD_MINNUMCELLS
    int FIELD_NPTISO_PISTON
    int FIELD_NPTISO_PDIFF
    int FIELD_PERIODIC
    int FIELD_SIMTIME
    int FIELD_MIN_GLOBAL_CUT
    int FIELD_THERMO_SWITCH
    int FIELD_THERMO_VIRTUAL
    int FIELD_TEMPERATURE
    int FIELD_LANGEVIN_GAMMA
    int FIELD_SWIMMING_PARTICLES_EXIST
    IF ROTATION:
        int FIELD_LANGEVIN_GAMMA_ROTATION
    IF NPT:
        int FIELD_NPTISO_G0
        int FIELD_NPTISO_GV
    int FIELD_MAX_OIF_OBJECTS

    void mpi_bcast_parameter(int p)

cdef extern from "communication.hpp":
    extern int n_nodes
    void mpi_set_time_step(double time_step)

cdef extern from "integrate.hpp":
    double time_step
    extern int integ_switch
    extern double sim_time
    extern double verlet_reuse
    extern double skin
    extern bool set_py_interrupt

cdef extern from "domain_decomposition.hpp":
    ctypedef struct  DomainDecomposition:
        int cell_grid[3]
        double cell_size[3]
        bool fully_connected[3]

    extern DomainDecomposition dd
    extern int max_num_cells
    extern int min_num_cells
    extern double max_skin
    int calc_processor_min_num_cells(const Vector3i & grid)


cdef extern from "particle_data.hpp":
    extern int n_part
    extern bool swimming_particles_exist

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    double dpd_gamma
    double dpd_r_cut
    extern double max_cut
    extern int max_seen_particle
    extern int max_seen_particle_type
    extern double max_cut_nonbonded
    extern double max_cut_bonded
    extern double min_global_cut


cdef extern from "thermostat.hpp":
    extern double nptiso_gamma0
    extern double nptiso_gammav
    extern double temperature
    extern int thermo_switch

cdef extern from "dpd.hpp":
    extern int dpd_wf
    extern double dpd_tgamma
    extern double dpd_tr_cut
    extern int dpd_twf


cdef extern from "cells.hpp":
    extern double max_range
    ctypedef struct CellStructure:
        int type
        bool use_verlet_list

    CellStructure cell_structure

cdef extern from "layered.hpp":
    extern int n_layers

cdef extern from "rattle.hpp":
    extern int n_rigidbonds


cdef extern from "tuning.hpp":
    extern int timing_samples


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

cdef extern from "immersed_boundaries.hpp":
    extern ImmersedBoundaries immersed_boundaries

cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int max_oif_objects

cdef extern from "forcecap.hpp":
    double forcecap_get()
    void forcecap_set(double forcecap)
