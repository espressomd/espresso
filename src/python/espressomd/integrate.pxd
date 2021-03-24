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
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map

include "myconfig.pxi"

cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp" nogil:
    cdef int python_integrate(int n_steps, cbool recalc_forces, int reuse_forces)
    cdef int mpi_steepest_descent(int max_steps)
    cdef void integrate_set_sd() except +
    cdef void integrate_set_nvt()
    cdef void integrate_set_steepest_descent(const double f_max, const double gamma,
                                             const double max_displacement) except +
    cdef extern cbool skin_set
    cdef extern cbool set_py_interrupt
    cdef void integrate_set_bd()

IF NPT:
    cdef extern from "integrate.hpp" nogil:
        cdef void integrate_set_npt_isotropic(double ext_pressure, double piston,
                                              cbool xdir_rescale, cbool ydir_rescale,
                                              cbool zdir_rescale, cbool cubic_box) except +

IF STOKESIAN_DYNAMICS:
    cdef extern from "stokesian_dynamics/sd_interface.hpp":
        void set_sd_viscosity(double eta) except +
        void set_sd_radius_dict(const unordered_map[int, double] & radius_dict) except +
        void set_sd_flags(int flg)

    cpdef enum flags:
        NONE = 0,
        SELF_MOBILITY = 1 << 0,
        PAIR_MOBILITY = 1 << 1,
        LUBRICATION = 1 << 2,
        FTS = 1 << 3

cdef inline int _integrate(int nSteps, cbool recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)
