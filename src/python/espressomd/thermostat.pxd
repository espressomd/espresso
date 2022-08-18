#
# Copyright (C) 2013-2022 The ESPResSo project
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
    int THERMO_SD

    cdef cppclass BaseThermostat:
        stdint.uint32_t rng_seed()
        stdint.uint64_t rng_counter()
        cbool is_seed_required()

    IF PARTICLE_ANISOTROPY:
        cdef cppclass LangevinThermostat(BaseThermostat):
            Vector3d gamma_rotation
            Vector3d gamma
        cdef cppclass BrownianThermostat(BaseThermostat):
            Vector3d gamma_rotation
            Vector3d gamma
    ELSE:
        cdef cppclass LangevinThermostat(BaseThermostat):
            double gamma_rotation
            double gamma
        cdef cppclass BrownianThermostat(BaseThermostat):
            double gamma_rotation
            double gamma
    cdef cppclass IsotropicNptThermostat(BaseThermostat):
        double gamma0
        double gammav
    cdef cppclass ThermalizedBondThermostat(BaseThermostat):
        pass
    IF DPD:
        cdef cppclass DPDThermostat(BaseThermostat):
            pass
    IF STOKESIAN_DYNAMICS:
        cdef cppclass StokesianThermostat(BaseThermostat):
            pass

    LangevinThermostat langevin
    BrownianThermostat brownian
    IsotropicNptThermostat npt_iso
    ThermalizedBondThermostat thermalized_bond
    IF DPD:
        DPDThermostat dpd
    IF STOKESIAN_DYNAMICS:
        StokesianThermostat stokesian

    void langevin_set_rng_seed(stdint.uint32_t seed)
    void brownian_set_rng_seed(stdint.uint32_t seed)
    void npt_iso_set_rng_seed(stdint.uint32_t seed)
    IF DPD:
        void dpd_set_rng_seed(stdint.uint32_t seed)
    IF STOKESIAN_DYNAMICS:
        void stokesian_set_rng_seed(stdint.uint32_t seed)

    void langevin_set_rng_counter(stdint.uint64_t counter)
    void brownian_set_rng_counter(stdint.uint64_t counter)
    void npt_iso_set_rng_counter(stdint.uint64_t counter)
    IF DPD:
        void dpd_set_rng_counter(stdint.uint64_t counter)
    IF STOKESIAN_DYNAMICS:
        void stokesian_set_rng_counter(stdint.uint64_t counter)

    IF PARTICLE_ANISOTROPY:
        void mpi_set_brownian_gamma(const Vector3d & gamma)
        void mpi_set_brownian_gamma_rot(const Vector3d & gamma)

        void mpi_set_langevin_gamma(const Vector3d & gamma)
        void mpi_set_langevin_gamma_rot(const Vector3d & gamma)
    ELSE:
        void mpi_set_brownian_gamma(const double & gamma)
        void mpi_set_brownian_gamma_rot(const double & gamma)

        void mpi_set_langevin_gamma(const double & gamma)
        void mpi_set_langevin_gamma_rot(const double & gamma)

    void mpi_set_thermo_virtual(cbool thermo_virtual)
    void mpi_set_temperature(double temperature)
    void mpi_set_thermo_switch(int thermo_switch)

    IF NPT:
        void mpi_set_nptiso_gammas(double gamma0, double gammav)

cdef extern from "stokesian_dynamics/sd_interface.hpp":
    IF STOKESIAN_DYNAMICS:
        void set_sd_kT(double kT) except +
        double get_sd_kT()

cdef extern from "grid_based_algorithms/lb_interface.hpp":
    double lb_lbfluid_get_kT "LB::get_kT"() except +

cdef extern from "grid_based_algorithms/lb_particle_coupling.hpp":
    void lb_lbcoupling_set_rng_state(stdint.uint64_t) except +
    stdint.uint64_t lb_lbcoupling_get_rng_state() except +
    void lb_lbcoupling_set_gamma(double) except +
    double lb_lbcoupling_get_gamma() except +
    cbool lb_lbcoupling_is_seed_required() except +
    void mpi_bcast_lb_particle_coupling()
