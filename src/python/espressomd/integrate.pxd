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
from __future__ import print_function, absolute_import
from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp cimport bool as cbool

cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp" nogil:
    cdef int python_integrate(int n_steps, int recalc_forces, int reuse_forces)
    cdef void integrate_set_nvt()
    cdef int integrate_set_npt_isotropic(double ext_pressure, double piston, int xdir, int ydir, int zdir, int cubic_box)
    cdef extern cbool skin_set
cdef inline int _integrate(int nSteps, int recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling::RuntimeError":
    cdef cppclass ErrorLevel:
        pass

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling::RuntimeError::ErrorLevel":
    cdef ErrorLevel WARNING
    cdef ErrorLevel ERROR

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling":
    cdef cppclass RuntimeError:
        string format()
        ErrorLevel level()

cdef extern from "errorhandling.hpp" namespace "ErrorHandling":
    cdef vector[RuntimeError]mpi_gather_runtime_errors()

cdef extern from "minimize_energy.hpp":
    void minimize_energy_init(const double f_max, const double gamma, const int max_steps, const double max_displacement)
cdef extern from "communication.hpp":
    int mpi_minimize_energy()


# cdef class Integrator:
#     cdef public _method
#     cdef public _steepest_descent_params
