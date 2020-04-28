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
from libcpp cimport bool as cbool
from libc cimport stdint

include "myconfig.pxi"
from .utils cimport Vector3d

cdef extern from "thermostat.hpp":
    double temperature
    int thermo_switch
    cbool thermo_virtual
    int THERMO_OFF
    int THERMO_LANGEVIN
    int THERMO_LB
    int THERMO_NPT_ISO
    int THERMO_DPD
    int THERMO_BROWNIAN

    IF PARTICLE_ANISOTROPY:
        ctypedef struct langevin_thermostat_struct "LangevinThermostat":
            Vector3d gamma_rotation
            Vector3d gamma
        ctypedef struct brownian_thermostat_struct "BrownianThermostat":
            Vector3d gamma_rotation
            Vector3d gamma
    ELSE:
        ctypedef struct langevin_thermostat_struct "LangevinThermostat":
            double gamma_rotation
            double gamma
        ctypedef struct brownian_thermostat_struct "BrownianThermostat":
            double gamma_rotation
            double gamma
    ctypedef struct npt_iso_thermostat_struct "IsotropicNptThermostat":
        double gamma0
        double gammav

    void langevin_set_rng_state(stdint.uint64_t counter)
    void brownian_set_rng_state(stdint.uint64_t counter)
    void npt_iso_set_rng_state(stdint.uint64_t counter)
    IF DPD:
        void dpd_set_rng_state(stdint.uint64_t counter)

    cbool langevin_is_seed_required()
    cbool brownian_is_seed_required()
    cbool npt_iso_is_seed_required()
    IF DPD:
        cbool dpd_is_seed_required()

    stdint.uint64_t langevin_get_rng_state()
    stdint.uint64_t brownian_get_rng_state()
    stdint.uint64_t npt_iso_get_rng_state()
    IF DPD:
        stdint.uint64_t dpd_get_rng_state()

cdef extern from "script_interface/Globals.hpp":
    # links intern C-struct with python object
    cdef extern langevin_thermostat_struct langevin
    cdef extern brownian_thermostat_struct brownian
    cdef extern npt_iso_thermostat_struct npt_iso

cdef extern from "npt.hpp":
    ctypedef struct nptiso_struct:
        double p_ext
        double p_inst
        double p_diff
        double piston
    extern nptiso_struct nptiso
