# Copyright (C) 2010-2018 The ESPResSo project
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
import unittest as ut
import numpy as np

import espressomd
from espressomd import has_features
if(has_features(["SWIMMER_REACTIONS"])):
    from espressomd.swimmer_reaction import Reaction


@ut.skipIf(not has_features(["SWIMMER_REACTIONS"]),
           "Features missing")
class ReactionTest(ut.TestCase):

    def test_reaction(self):
        system = espressomd.System(box_l=[10.0, 10.0, 10.0])
        system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

        system.cell_system.skin = 0.1
        system.time_step = 0.01

        system.part.add(pos=0.5 * system.box_l, type=0)
        reactant = system.part.add(
            pos=0.5 * system.box_l - [0, 0, 1], type=1, v=[0, 0, 0])
        product = system.part.add(
            pos=0.5 * system.box_l + [0, 0, 1], type=2, v=[0, 0, 0])

        # Set up the reaction, rate is infinity to almost surely trigger a
        # reaction
        Reaction(reactant_type=1, catalyzer_type=0,
                 product_type=2, ct_range=3.0, ct_rate=np.inf)

        system.integrator.run(
            2)  # run two reactions and MD. this also tests that the particles are not switched back after the first reaction (which should switch the particles)

        # Check that particles have switched type (due to the first reaction)
        # the second reaction may not change the particles back since they are
        # now not in the right orientation for a reaction
        self.assertEqual(reactant.type, 2)
        self.assertEqual(product.type, 1)

if __name__ == '__main__':
    ut.main()
