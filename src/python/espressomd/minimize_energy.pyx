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


def steepest_descent(system, *args, **kwargs):
    """
    Steepest descent algorithm for energy minimization. Wrapper for
    :class:`espressomd.integrate.SteepestDescent`.

    Parameters
    ----------
    system : :obj:`espressomd.system.System`
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

    Returns
    -------
    :obj:`bool`
        Whether the steepest descent has converged.

    """
    cdef object old_integrator
    old_integrator = system.integrator.get_state()
    system.integrator.set_steepest_descent(*args, **kwargs)
    steps = system.integrator.run(kwargs['max_steps'])
    system.integrator.__setstate__(old_integrator)
    return steps < kwargs['max_steps']
