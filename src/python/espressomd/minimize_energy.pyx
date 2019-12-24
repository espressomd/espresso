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

from espressomd.utils import is_valid_type

cdef class MinimizeEnergy:
    """
    Steepest descent algorithm for energy minimization. Wrapper for
    :class:`espressomd.integrate.SteepestDescent`.

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
    cdef object _integrator_handle
    cdef object _old_integrator
    cdef int _max_steps

    def __getstate__(self):
        return {'integrator_handle': self._integrator_handle,
                'old_integrator': self._old_integrator,
                'max_steps': self._max_steps}

    def __setstate__(self, state):
        self._integrator_handle = state['integrator_handle']
        self._old_integrator = state['old_integrator']
        self._max_steps = state['max_steps']

    def __init__(self, integrator_handle):
        self._integrator_handle = integrator_handle

    def init(self, *args, **kwargs):
        """
        Initialize the steepest descent integrator.

        """
        self._old_integrator = self._integrator_handle.get_state()
        self._max_steps = kwargs['max_steps']
        self._integrator_handle.set_steepest_descent(*args, **kwargs)

    def minimize(self):
        """
        Perform energy minimization sweep.

        """
        self._integrator_handle.run(self._max_steps)

    def disable(self, *args, **kwargs):
        """
        Restore the original integrator.

        """
        self._integrator_handle.__setstate__(self._old_integrator)
