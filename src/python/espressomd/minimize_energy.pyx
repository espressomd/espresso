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

from . cimport minimize_energy
from espressomd.utils import is_valid_type

cdef class MinimizeEnergy:
    """
    Steepest descent algorithm for energy minimization.

    Particles located at :math:`\\vec{r}_i` at integration step :math:`i` and
    experiencing a potential :math:`\mathcal{H}(\\vec{r}_i)` are displaced
    according to the equation:

    :math:`\\vec{r}_{i+1} = \\vec{r}_i - \\gamma\\nabla\mathcal{H}(\\vec{r}_i)`

    Parameters
    ----------
    f_max : :obj:`float`
        Convergence criterion. Minimization stops when the maximal force on
        particles in the system is lower than this threshold. Set this to 0
        when running minimization in a loop that stops when a custom
        convergence criterion is met.
    gamma : :obj:`float`
        Dampening constant.
    max_steps : :obj:`int`
        Maximal number of iterations.
    max_displacement : :obj:`float`
        Maximal allowed displacement per step. Typical values for a LJ liquid
        are in the range of 0.1% to 10% of the particle sigma.

    """
    cdef object _params

    def __getstate__(self):
        return self._params

    def __setstate__(self, params):
        self._params = params

    def __init__(self, *args, **kwargs):
        if len(args) == 0:
            # Initialize default values
            self._params = self.default_params()
            return

        # Check if all required keys are given
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        self._params = kwargs
        self.validate_params()

    def init(self, *args, **kwargs):
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        self._params = kwargs
        self.validate_params()

    def default_params(self):
        para = dict()
        para["f_max"] = 0.0
        para["gamma"] = 0.0
        para["max_steps"] = 0
        para["max_displacement"] = 0.0
        return para

    def required_keys(self):
        return "f_max", "gamma", "max_steps", "max_displacement"

    def validate_params(self):
        if self._params["f_max"] < 0:
            raise ValueError(
                "f_max has to be a positive floating point number")
        if self._params["gamma"] < 0:
            raise ValueError(
                "gamma has to be a positive floating point number")
        if self._params["max_steps"] < 0 or not is_valid_type(
                self._params["max_steps"], int):
            raise ValueError(
                "max_steps has to be a positive integer")
        if self._params["max_displacement"] < 0:
            raise ValueError(
                "max_displacement has to be a positive floating point number")

    def minimize(self):
        """
        Perform energy minimization sweep.

        """
        minimize_energy_init(self._params["f_max"], self._params["gamma"],
                             self._params["max_steps"],
                             self._params["max_displacement"])
        mpi_minimize_energy()
