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
from libcpp.string cimport string  # import std::string as string
from libcpp.list cimport list  # import std::list as list
from libcpp cimport bool as cbool


cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp":
    cdef int python_integrate(int n_steps, int recalc_forces, int reuse_forces)
    cdef void integrate_set_nvt()
    cdef int integrate_set_npt_isotropic(double ext_pressure, double piston, int xdir, int ydir, int zdir, int cubic_box)
    cdef extern cbool skin_set

cdef extern from "errorhandling.hpp":
    cdef list[string] mpiRuntimeErrorCollectorGather()
