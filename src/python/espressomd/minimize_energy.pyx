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
# Minimize Energy

from __future__ import print_function, absolute_import
from . cimport minimize_energy
from espressomd.utils import is_valid_type

cdef class MinimizeEnergy(object):
    """
    Initialize steepest descent energy minimization.

    Parameters
    ----------
    f_max : :obj:`float`
            Maximal allowed force.
    gamma : :obj:`float`
            Dampening constant.
    max_steps : :obj:`int`
                Maximal number of iterations.
    max_displacement : :obj:`float`
                       Maximal allowed displacement per step.

    """
    cdef object _params

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
        if self._params["max_steps"] < 0 or not is_valid_type(self._params["max_steps"], int):
            raise ValueError(
                "max_steps has to be a positive integer")
        if self._params["max_displacement"] < 0:
            raise ValueError(
                "max_displacement has to be a positive floating point number")

    def minimize(self):
        """
        Perform energy minimization sweep.

        """
        minimize_energy_init(self._params["f_max"], self._params["gamma"], self._params[
                             "max_steps"], self._params["max_displacement"])
        mpi_minimize_energy()
