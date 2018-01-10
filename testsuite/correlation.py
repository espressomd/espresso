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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from numpy.random import random
from espressomd.interactions import FeneBond
from espressomd.observables import *
from espressomd.correlators import *
import pickle

class Observables(ut.TestCase):
    # Handle for espresso system
    es = espressomd.System()

    def test_corr(self):
        s = self.es
        s.box_l = 10, 10, 10
        s.cell_system.skin = 0.4
        # s.periodicity=0,0,0
        s.time_step = 0.01
        s.thermostat.turn_off()
        s.part.add(id=0, pos=(0, 0, 0), v=(1, 2, 3))

        O = ParticlePositions(ids=(0,))
        C2 = Correlator(obs1=O, dt=0.01, tau_lin=10, tau_max=10.0,
                        corr_operation="square_distance_componentwise")

        s.integrator.run(1000)
        s.auto_update_correlators.add(C2)
        s.integrator.run(20000)


        corr = C2.result()

        # Check pickling
        C_unpickeled = pickle.loads(pickle.dumps(C2))
        self.assertTrue(np.array_equal(corr, C_unpickeled.result()))

        for i in range(corr.shape[0]):
            t = corr[i, 0]
            self.assertAlmostEqual(corr[i, 2], t * t, places=3)
            self.assertAlmostEqual(corr[i, 3], 4 * t * t, places=3)
            self.assertAlmostEqual(corr[i, 4], 9 * t * t, places=3)


if __name__ == "__main__":
    ut.main()
