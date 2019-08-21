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
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class MMM2D_setup(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.01

    def test_nonzero_maxcut(self):

        self.system.cell_system.set_layered(
            n_layers=20, use_verlet_lists=False)

        self.system.periodicity = [1, 1, 0]

        self.system.part.add(pos=(1.0, 1.0, 1.0), q=1)
        self.system.part.add(pos=(9.0, 9.0, 9.0), q=-1)

        # No interaction so far: max_cut == 0
        self.assertTrue(self.system.cell_system.get_state()["max_cut"] == 0)

        mmm2d = espressomd.electrostatics.MMM2D(prefactor=1.0, maxPWerror=1e-3)
        self.system.actors.add(mmm2d)
        
        # MMM2D alters max_cut to small eps
        self.assertTrue(self.system.cell_system.get_state()["max_cut"] != 0)
       
        # Test integration
        self.system.integrator.run(0)

if __name__ == "__main__":
    ut.main()
