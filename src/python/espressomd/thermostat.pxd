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

cdef extern from "Vector.hpp":
    cppclass Vector3d:
        double & operator[](int i)

cdef extern from "thermostat.hpp":
    double temperature
    int thermo_switch
    int THERMO_OFF
    int THERMO_LANGEVIN
    int THERMO_LB
    int THERMO_NPT_ISO
    int THERMO_DPD
    int THERMO_INTER_DPD

    IF PARTICLE_ANISOTROPY:
        Vector3d langevin_gamma_rotation
        Vector3d langevin_gamma
    ELSE:
        double langevin_gamma_rotation
        double langevin_gamma
