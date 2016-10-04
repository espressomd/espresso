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
include "myconfig.pxi"

cdef extern from "communication.hpp":
    void mpi_bcast_parameter(int p)

cdef extern from "global.hpp":
    int FIELD_TEMPERATURE
    int FIELD_THERMO_SWITCH
    int FIELD_TEMPERATURE
    int FIELD_LANGEVIN_GAMMA
    IF ROTATION:
        int FIELD_LANGEVIN_GAMMA_ROTATION
    IF NPT:
        int FIELD_NPTISO_G0
        int FIELD_NPTISO_GV

cdef extern from "thermostat.hpp":
    double temperature
    int thermo_switch
    double langevin_gamma
    int THERMO_OFF
    int THERMO_LANGEVIN
    IF ROTATION:
        IF ROTATIONAL_INERTIA:
            double langevin_gamma_rotation[3]
        ELSE:
            double langevin_gamma_rotation
    int THERMO_LB
    int THERMO_NPT_ISO
    int THERMO_DPD
    int THERMO_INTER_DPD
    IF ROTATIONAL_INERTIA:
        double langevin_gamma_rotation[3]
    ELSE:
        double langevin_gamma_rotation

