#
# Copyright (C) 2013-2020 The ESPResSo project
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
import sys
import numpy as np
import unittest as ut
import unittest_decorators as utx
from espressomd import system, minimize_energy, generic_dd
from espressomd.interactions import FeneBond


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class Generic_DD_Metric(ut.TestCase):
    """"Test metric functionality of generic_dd."""
    box = 16.0
    s = system.System(box_l=[box, box, box])
    nproc = s.cell_system.get_state()["n_nodes"]
    F = FeneBond(k=1, d_r_max=2)
    s.bonded_inter.add(F)

    @classmethod
    def setUpClass(cls):
        cls.s.time_step = .000001
        cls.s.cell_system.skin = 0.0
        cls.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.5, shift="auto")
        cls.dd = cls.s.cell_system.set_generic_dd("cart")

    def test_linear_combination(self):
        m = self.dd.metric("ncells")
        m2 = self.dd.metric("-ncells+2.0*ncells")
        self.assertEqual(m.average(), m2.average())
        self.assertEqual(m.maximum(), m2.maximum())
        self.assertEqual(m.imbalance(), m2.imbalance())

    def test_zero_without_particles(self):
        self.s.part.clear()
        m = self.dd.metric("npart")
        self.assertEqual(m.average(), 0.0)
        self.assertEqual(m.maximum(), 0.0)
        m = self.dd.metric("ndistpairs")
        self.assertEqual(m.average(), 0.0)
        self.assertEqual(m.maximum(), 0.0)
        m = self.dd.metric("nbondedia")
        self.assertEqual(m.average(), 0.0)
        self.assertEqual(m.maximum(), 0.0)
        m = self.dd.metric("npart+2.0*ndistpairs+3.0*nbondedia")
        self.assertEqual(m.average(), 0.0)
        self.assertEqual(m.maximum(), 0.0)

    def test_npart(self):
        self.s.part.clear()
        self.s.part.add(pos=[1,1,1])
        m = self.dd.metric("npart")
        self.assertEqual(m.average(), 1.0 / self.nproc)
        self.assertEqual(m.maximum(), 1.0)
        self.assertEqual(m.imbalance(), self.nproc)

    def test_ndistpairs(self):
        self.s.part.clear()
        self.s.part.add(pos=[1,1,1])
        m = self.dd.metric("ndistpairs")
        self.assertEqual(m.average(), 0.0)
        self.assertEqual(m.maximum(), 0.0)
        self.s.part.add(pos=[1,1,2])
        self.assertEqual(m.average(), 1.0 / self.nproc)
        self.assertEqual(m.maximum(), 1.0)
        self.assertEqual(m.imbalance(), self.nproc)
        self.s.part.add(pos=[1,3,1])
        self.assertEqual(m.average(), 3.0 / self.nproc)
        self.s.part.add(pos= 16 * np.random.rand(100,3))
        self.s.integrator.run(1) # Force a ghost exchange
        ground_truth = len(self.s.cell_system.get_pairs_(np.max(self.s.box_l)))
        self.assertEqual(m.average(),  ground_truth / self.nproc)

    def test_nbondedia(self):
        self.s.part.clear()
        self.s.part.add(pos=[1,1,1])
        self.s.part.add(pos=[2,2,2])
        self.s.part.add(pos=[3,3,3])
        self.s.part[0].add_bond((self.F, 1))
        self.s.part[0].add_bond((self.F, 2))
        m = self.dd.metric("nbondedia")
        self.assertEqual(m.average(), 2.0 / self.nproc)
        self.assertEqual(m.maximum(), 2.0) # Both bonds are stored on particle 0


if __name__ == "__main__":
    ut.main()
