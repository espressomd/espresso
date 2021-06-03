#
# Copyright (C) 2013-2021 The ESPResSo project
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
This script uses the System class from the :file:`gibbs_ensemble.py`
script and creates a subclass in which you can create your system, add
energy corrections etc.
"""

import espressomd
import gibbs_ensemble
import numpy as np

espressomd.assert_features("LENNARD_JONES")

# Lennard Jones parameters
LJ_EPSILON = 1.0
LJ_SIGMA = 1.0
LJ_CUTOFF = 2.5 * LJ_SIGMA
LJ_SHIFT = - (np.power(LJ_SIGMA / LJ_CUTOFF, 12) -
              np.power(LJ_SIGMA / LJ_CUTOFF, 6))


def calc_shift_correction(density):
    """ Shift correction by integrating the constant part from 0 to the cutoff
        radius """
    return 4. / 3. * np.pi * LJ_CUTOFF**3 * LJ_SHIFT * 4 * LJ_EPSILON * density


def calc_tail_correction(density):
    """ Lennard-Jones truncation correction (See Frenkel, Smit, Understanding
        molecular simulation, eq. (3.2.5)) """
    return 8. / 3. * np.pi * density * LJ_EPSILON * LJ_SIGMA**3 * \
        (1. / 3. * np.power(LJ_SIGMA / LJ_EPSILON, 9) -
         np.power(LJ_SIGMA / LJ_EPSILON, 3))


class Gibbs_Client(gibbs_ensemble.Client):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def energy(self):
        """ Additional correction terms """
        return super().energy  # + \
        # calc_tail_correction(self.density) + \
        # calc_shift_correction(self.density)

    def init_system(self):
        """ Initialize the system (this is executed only when the run
            method of the process is executed) """
        self.system = set_up_system(self.init_box_length, self.init_num_part)


def set_up_system(box_length, num_part):
    """ Use this function to create the system """

    system = espressomd.System(box_l=3 * [box_length])

    system.cell_system.set_n_square(use_verlet_lists=False)
    system.cell_system.skin = 0
    system.time_step = 0.01

    system.part.add(pos=np.random.random((num_part, 3)) * box_length)

    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=LJ_EPSILON,
        sigma=LJ_SIGMA,
        cutoff=LJ_CUTOFF,
        shift=LJ_SHIFT)

    return system
