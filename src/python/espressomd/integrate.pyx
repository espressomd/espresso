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
include "myconfig.pxi"
from utils cimport *

cdef int _integrate(int nSteps, int recalc_forces, int reuse_forces):
    with nogil:
        return python_integrate(nSteps, recalc_forces, reuse_forces)

def handle_errors(msg):
    errors = mpi_gather_runtime_errors()
    for err in errors:
        print(err.format())

    for err in errors:
    # Cast because cython does not support typed enums completely
        if <int> err.level() == <int> ERROR:
            raise Exception(msg)
    

def integrate(nSteps, recalc_forces=False, reuse_forces=False):
    """integrate(nSteps, recalc_forces=False, reuse_forces=False)"""

    check_type_or_throw_except(
        nSteps, 1, int, "Integrate requires a positive integer for the number of steps")
    check_type_or_throw_except(
        recalc_forces, 1, bool, "recalc_forces has to be a bool")
    check_type_or_throw_except(
        reuse_forces, 1, bool, "reuse_forces has to be a bool")
    
    if (_integrate(nSteps, recalc_forces, reuse_forces)):
        handle_errors("Encoutered errors during integrate")

def set_integrator_nvt():
    integrate_set_nvt()


def set_integrator_isotropic_npt(ext_pressure=0.0, piston=0.0, xdir=0, ydir=0, zdir=0, cubic_box=False):
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
        handle_errors("Encoutered errors setting up the NPT integrator")
