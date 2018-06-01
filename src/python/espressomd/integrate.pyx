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
import espressomd.code_info
from espressomd.utils cimport *

cdef class Integrator(object):
    """
    Integrator class.

    This class interfaces the Velocity Verlet integrator.

    """

    cdef str _method
    cdef object _steepest_descent_params

    def __init__(self):
        self._method = "VV"
        self._steepest_descent_params = {}

    def run(self, steps=1, recalc_forces=False, reuse_forces=False):
        """
        Run the integrator.

        Parameters
        ----------
        steps : :obj:`int`
            Number of time steps to integrate.
        recalc_forces : :obj:`bool`, optional
            Recalculate the forces regardless of whether they are reusable.
        reuse_forces : :obj:`bool`, optional
            Reuse the forces from previous time step.

        """

        if self._method == "VV":
            check_type_or_throw_except(
                steps, 1, int, "Integrate requires a positive integer for the number of steps")
            check_type_or_throw_except(
                recalc_forces, 1, bool, "recalc_forces has to be a bool")
            check_type_or_throw_except(
                reuse_forces, 1, bool, "reuse_forces has to be a bool")

            if (_integrate(steps, recalc_forces, reuse_forces)):
                handle_errors("Encoutered errors during integrate")
        elif self._method == "STEEPEST_DESCENT":
            minimize_energy_init(self._steepest_descent_params["f_max"],
                                 self._steepest_descent_params["gamma"],
                                 steps,
                                 self._steepest_descent_params["max_displacement"])
            mpi_minimize_energy()

    def set_steepest_descent(self, *args, **kwargs):
        """
        Set parameters for steepest descent.

        .. seealso::
            :class:`espressomd.minimize_energy.MinimizeEnergy`

        """

        req = ["f_max", "gamma", "max_displacement"]
        for key in kwargs:
            if not key in req:
                raise Exception("Set required parameter %s first." %key)

        self._steepest_descent_params.update(kwargs)
        self._method = "STEEPEST_DESCENT"

    def set_vv(self):
        """
        Set the integration method to Velocity Verlet.

        """
        self._method = "VV"

    def set_nvt(self):
        """
        Set the integration method to NVT.

        """

        integrate_set_nvt()

    def set_isotropic_npt(self, ext_pressure, piston, direction=[0,0,0], cubic_box=False):
        """
        Set the integration method to NPT.

        Parameters
        ----------
        ext_pressure : :obj:`float`
            The external pressure.
        piston : :obj:`float`
            The mass of the applied piston.
        direction : :obj:`list`, optional
            Three integers to set the box geometry for non-cubic boxes
        cubic_box : :obj:`bool`, optional
            If this optional parameter is true, a cubic box is assumed.

        """

        if "NPT" not in espressomd.code_info.features():
            raise Exception("NPT is not compiled in")
        check_type_or_throw_except(
            ext_pressure, 1, float, "NPT parameter ext_pressure must be a float")
        check_type_or_throw_except(
            piston, 1, float, "NPT parameter piston must be a float")
        check_type_or_throw_except(
            direction, 3, int, "NPT parameter direction must be an three ints")
        if (integrate_set_npt_isotropic(ext_pressure, piston, direction[0], direction[1], direction[2], cubic_box)):
            handle_errors("Encoutered errors setting up the NPT integrator")
