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
include "myconfig.pxi"
from libcpp cimport bool

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
    int FIELD_BROWNIAN_GAMMA
    IF ROTATION:
        int FIELD_LANGEVIN_GAMMA_ROTATION
        int FIELD_BROWNIAN_GAMMA_ROTATION
    IF NPT:
        int FIELD_NPTISO_G0
        int FIELD_NPTISO_GV
    int FIELD_MAX_OIF_OBJECTS

    void mpi_bcast_parameter(int p)

cdef extern from "communication.hpp":
    void mpi_set_time_step(double time_step) except +

cdef extern from "integrate.hpp":
    double time_step
    extern int integ_switch
    extern double sim_time
    extern double verlet_reuse
    extern double skin

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    extern int max_seen_particle_type
    extern double min_global_cut
    double maximal_cutoff_bonded()
    double maximal_cutoff_nonbonded()

cdef extern from "rattle.hpp":
    extern int n_rigidbonds

cdef extern from "tuning.hpp":
    extern int timing_samples

cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int max_oif_objects

cdef extern from "forcecap.hpp":
    double forcecap_get()
    void forcecap_set(double forcecap)
