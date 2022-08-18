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

include "myconfig.pxi"

cdef extern from "integrate.hpp":
    double get_time_step()

cdef extern from "integrate.hpp" nogil:
    cdef int python_integrate(int n_steps, cbool recalc_forces, int reuse_forces)
    cdef int mpi_steepest_descent(int max_steps)
    cdef extern cbool set_py_interrupt

cdef inline int _integrate(int nSteps, cbool recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)
