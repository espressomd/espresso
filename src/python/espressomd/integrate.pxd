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
from libc cimport stdint
from libcpp cimport bool as cbool
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map

include "myconfig.pxi"

cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp" nogil:
    cdef int python_integrate(int n_steps, cbool recalc_forces, int reuse_forces)
    cdef void integrate_set_sd()
    cdef void integrate_set_nvt()
    cdef int integrate_set_steepest_descent(const double f_max, const double gamma,
                                            const int max_steps, const double max_displacement)
    cdef extern cbool skin_set
    cdef extern cbool set_py_interrupt
    cdef void integrate_set_bd()

    stdint.uint64_t get_integrator_counter()
    void set_integrator_counter(stdint.uint64_t value)

IF NPT:
    cdef extern from "integrate.hpp" nogil:
        cdef int integrate_set_npt_isotropic(double ext_pressure, double piston,
                                             cbool xdir_rescale, cbool ydir_rescale,
                                             cbool zdir_rescale, cbool cubic_box)

cdef extern from "stokesian_dynamics/sd_interface.hpp":
    IF(STOKESIAN_DYNAMICS or STOKESIAN_DYNAMICS_GPU):
        void set_sd_viscosity(double eta)
        double get_sd_viscosity()

        void set_sd_device(const string & dev)
        string get_sd_device()

        void set_sd_radius_dict(const unordered_map[int, double] & radius_dict)
        unordered_map[int, double] get_sd_radius_dict()

        void set_sd_flags(int flg)
        int get_sd_flags()

IF(STOKESIAN_DYNAMICS or STOKESIAN_DYNAMICS_GPU):
    cpdef enum flags:
        NONE = 0,
        SELF_MOBILITY = 1 << 0,
        PAIR_MOBILITY = 1 << 1,
        LUBRICATION = 1 << 2,
        FTS = 1 << 3

cdef inline int _integrate(int nSteps, cbool recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)

cdef extern from "communication.hpp":
    int mpi_steepest_descent(int max_steps)
