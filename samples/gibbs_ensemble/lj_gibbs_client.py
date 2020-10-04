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
"""
Client part of the Gibbs ensemble simulation. This script handles the
simulation boxes and communicates the energies to the host. The Monte-Carlo
part of the simulation is done by the :file:`gibbs_ensemble_socket.py` script.
"""

import argparse

import numpy as np

import espressomd
import gibbs

espressomd.assert_features("LENNARD_JONES")

seed = None

LJ_EPSILON = 1.0
LJ_SIGMA = 1.0
LJ_CUTOFF = 2.5
LJ_SHIFT = 0


parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int)
parser.add_argument('-H', '--host', type=str, default='localhost')
parser.add_argument('-p', '--port', type=int, required=True)
args = parser.parse_args()


if args.seed:
    np.random.seed(args.seed)


def calc_tail_correction(density, lj_epsilon, lj_sigma, lj_cutoff):
    '''
    Calculates the tail correction to the energies of the box.
    eq 3.2.5
    '''
    return 8.0 / 3.0 * np.pi * density * lj_epsilon * \
        lj_sigma**3 * (1.0 / 3.0 * np.power(lj_cutoff / lj_sigma, -9) -
                       np.power(lj_cutoff / lj_sigma, -3))


def calc_shift_correction(density, lj_epsilon, lj_cutoff, lj_shift):
    '''
    Calculates the shift correction to the energies of the box.
    difference in the potential integrated from 0 to cutoff distance
    '''
    return -8.0 / 3.0 * np.pi * density * \
        lj_epsilon * np.power(lj_cutoff, 3) * 4.0 * lj_shift


class LJGibbsClient(gibbs.Client):
    def __init__(self, system, lj_epsilon, lj_sigma, lj_cutoff, lj_shift):
        self._lj_epsilon = lj_epsilon
        self._lj_sigma = lj_sigma
        self._lj_cutoff = lj_cutoff
        self._lj_shift = lj_shift
        super().__init__(system)

    def density(self):
        return len(self._system.part) / self._system.volume()

    def energy(self):
        """Add Lennard-Jones shift and tail correction to the energy"""
        return super().energy() + \
            calc_tail_correction(self.density(), self._lj_epsilon, self._lj_sigma, self._lj_cutoff) + \
            calc_shift_correction(self.density(), self._lj_epsilon,
                                  self._lj_cutoff, self._lj_shift)



# init system
# The intial system has to be big enough to fit twice the
# maximum interaction range
box_l = 3 * LJ_CUTOFF
system = espressomd.System(box_l=[box_l, box_l, box_l])
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=LJ_EPSILON,
                                                       sigma=LJ_SIGMA,
                                                       cutoff=LJ_CUTOFF,
                                                       shift=LJ_SHIFT)
system.cell_system.set_n_square(use_verlet_lists=False)
system.cell_system.skin = 0

gibbs_client = LJGibbsClient(system, LJ_EPSILON, LJ_SIGMA, LJ_CUTOFF, LJ_SHIFT)
gibbs_client.run(args.host, args.port)
