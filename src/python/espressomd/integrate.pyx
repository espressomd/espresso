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
from espressomd.utils cimport *

cdef class Integrator:
    cdef str _method
    cdef object _steepest_descent_params

    def __init__(self):
        self._method = "VV"
        self._steepest_descent_params = {}

    def run(self, steps=1, recalc_forces=False, reuse_forces=False):
        if self._method == "VV":
            check_type_or_throw_except(
                steps, 1, int, "Integrate requires a positive integer for the number of steps")
            check_type_or_throw_except(
                recalc_forces, 1, bool, "recalc_forces has to be a bool")
            check_type_or_throw_except(
                reuse_forces, 1, bool, "reuse_forces has to be a bool")

            if (_integrate(steps, recalc_forces, reuse_forces)):
                self.handle_errors("Encoutered errors during integrate")
        elif self._method == "STEEPEST_DESCENT":
            minimize_energy_init(self._steepest_descent_params["f_max"], self._steepest_descent_params["gamma"], steps, self._steepest_descent_params["max_displacement"])
            minimize_energy()

    def set_steepest_descent(self, *args, **kwargs):
        req = ["f_max", "gamma", "max_displacement"]
        for key in kwargs:
            if not key in req:
                raise Exception("Set required parameter %s first." %key)
            
        self._steepest_descent_params.update(kwargs)
        self._method = "STEEPEST_DESCENT"

    cdef set_vv(self):
        self._method = "VV"

    cdef set_nvt(self):
        integrate_set_nvt()

    cdef set_isotropic_npt(self, ext_pressure=0.0, piston=0.0, xdir=0, ydir=0, zdir=0, cubic_box=False):
        IF NPT != 1:
            raise Exception("NPT is not compiled in")
        check_type_or_throw_except(
            ext_pressure, 1, float, "NPT parameter ext_pressure must be a float")
        check_type_or_throw_except(
            piston, 1, float, "NPT parameter piston must be a float")
        check_type_or_throw_except(
            xdir, 1, int, "NPT parameter xdir must be an int")
        check_type_or_throw_except(
            ydir, 1, int, "NPT parameter ydir must be an int")
        check_type_or_throw_except(
            zdir, 1, int, "NPT parameter zdir must be an int")
        if (integrate_set_npt_isotropic(ext_pressure, piston, xdir, ydir, zdir, cubic_box)):
            self.handle_errors("Encoutered errors setting up the NPT integrator")

            
    cdef handle_errors(self, msg):
        errors = mpi_gather_runtime_errors()
        for err in errors:
            print(err.format())

        for err in errors:
        # Cast because cython does not support typed enums completely
            if <int> err.level() == <int> ERROR:
                raise Exception(msg)




