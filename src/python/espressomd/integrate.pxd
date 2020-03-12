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

include "myconfig.pxi"

cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp" nogil:
    cdef int python_integrate(int n_steps, cbool recalc_forces, int reuse_forces)
    cdef void integrate_set_nvt()
    cdef int integrate_set_steepest_descent(const double f_max, const double gamma,
                                            const int max_steps, const double max_displacement)
    cdef extern cbool skin_set
    cdef extern cbool set_py_interrupt
    cdef void integrate_set_bd()

IF NPT:
    cdef extern from "integrate.hpp" nogil:
        cdef int integrate_set_npt_isotropic(double ext_pressure, double piston,
                                             cbool xdir_rescale, cbool ydir_rescale,
                                             cbool zdir_rescale, cbool cubic_box)

cdef inline int _integrate(int nSteps, cbool recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)

cdef extern from "communication.hpp":
    int mpi_steepest_descent(int max_steps)
