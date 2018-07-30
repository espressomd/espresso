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
# Integration test for exclusions
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np


@ut.skipIf(not espressomd.has_features(['EXCLUSIONS']), "Skipping test")
class Exclusions(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.s.part.clear()
        self.s.box_l = 3 * [10]
        self.s.cell_system.skin = 0.4
        self.s.time_step = 0.01

    def test_add_remove(self):
        self.s.part.add(id=0, pos=[0, 0, 0])
        self.s.part.add(id=1, pos=[0, 0, 0])
        self.s.part.add(id=2, pos=[0, 0, 0])

        self.s.part[0].add_exclusion(1)
        self.s.part[0].add_exclusion(2)
        self.assertTrue( (self.s.part[0].exclusions == [1, 2]).all() )
        self.s.part[0].delete_exclusion(1)
        self.assertEqual(self.s.part[0].exclusions, [2])
        self.s.part[0].delete_exclusion(2)
        self.assertEqual( list(self.s.part[0].exclusions), [])

    def test_transfer(self):
        self.s.part.add(id=0, pos=[0, 0, 0], v=[1., 1., 1])
        self.s.part.add(id=1, pos=[0, 0, 0])
        self.s.part.add(id=2, pos=[0, 0, 0])
        self.s.part.add(id=3, pos=[0, 0, 0])

        self.s.part[0].exclusions = [1, 2, 3]

        for i in range(15):
            self.s.integrator.run(100)
            self.assertTrue( (self.s.part[0].exclusions == [1, 2, 3]).all() )

    @ut.skipIf(not espressomd.has_features(['LENNARD_JONES']), "Skipping test")
    def test_particle_property(self):
        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=2.,
            cutoff=1.5, shift=0.0)

        self.s.part.add(id=0, pos=[0, 0, 0], type=0)
        self.s.part.add(id=1, pos=[1, 0, 0], type=0)

        pair_energy = self.s.analysis.energy()['total']
        self.assertGreater(pair_energy, 0.)

        pair_pressure = self.s.analysis.pressure()['total']
        self.assertGreater(pair_pressure, 0.)

        self.s.integrator.run(0)
        pair_force = self.s.part[0].f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(self.s.part[1].f[0], -pair_force, places=7)

        self.s.part.add(id=2, pos=[2, 0, 0], type=0)
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.analysis.energy()[
                               'total'], 2 * pair_energy)
        self.assertAlmostEqual(self.s.analysis.pressure()[
                               'total'], 2 * pair_pressure)
        self.assertAlmostEqual(self.s.part[2].f[0], -pair_force, places=7)

        self.s.part[1].exclusions = [0, 2]
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.analysis.energy()['total'], 0)
        self.assertAlmostEqual(self.s.analysis.pressure()['total'], 0)
        self.assertAlmostEqual(self.s.part[0].f[0], 0, places=7)
        self.assertAlmostEqual(self.s.part[1].f[0], 0, places=7)
        self.assertAlmostEqual(self.s.part[2].f[0], 0, places=7)

        self.s.part[1].exclusions = [0]
        self.assertAlmostEqual(self.s.analysis.energy()['total'], pair_energy)
        self.assertAlmostEqual(self.s.analysis.pressure()['total'], pair_pressure)
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.part[0].f[0], 0, places=7)
        self.assertAlmostEqual(self.s.part[1].f[0], pair_force, places=7)
        self.assertAlmostEqual(self.s.part[2].f[0], -pair_force, places=7)

        self.s.part[1].exclusions = []
        self.assertAlmostEqual(self.s.analysis.energy()[
                               'total'], 2 * pair_energy)
        self.assertAlmostEqual(self.s.analysis.pressure()[
                               'total'], 2 * pair_pressure)
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.part[0].f[0], pair_force, places=7)
        self.assertAlmostEqual(self.s.part[1].f[0], 0, places=7)
        self.assertAlmostEqual(self.s.part[2].f[0], -pair_force, places=7)

        self.s.part[1].exclusions = [0]
        self.assertAlmostEqual(self.s.analysis.energy()['total'], pair_energy)
        self.assertAlmostEqual(self.s.analysis.pressure()['total'], pair_pressure)
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.part[0].f[0], 0, places=7)
        self.assertAlmostEqual(self.s.part[1].f[0], pair_force, places=7)
        self.assertAlmostEqual(self.s.part[2].f[0], -pair_force, places=7)

    @ut.skipIf(not espressomd.has_features(['P3M']), "Skipping test")
    def test_electrostatics_not_excluded(self):
        from espressomd.electrostatics import P3M
        self.s.part.add(id=0, pos=[0, 0, 0], type=0, q=+1.)
        self.s.part.add(id=1, pos=[1, 0, 0], type=0, q=-1.)

        # Small alpha means large short-range contribution
        self.s.actors.add(P3M(prefactor=1, r_cut=3.0, accuracy=1e-3,
                                  mesh=32, cao=7, alpha=0.1, tune=False))

        # Only short-range part of the coulomb energy
        pair_energy = self.s.analysis.energy()[('coulomb', 0)]
        self.assertGreater(abs(pair_energy), 0.)

        self.s.integrator.run(0)
        pair_force = self.s.part[0].f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(self.s.part[1].f[0], -pair_force, places=7)

        pair_pressure = self.s.analysis.pressure()[('coulomb', 0)]
        self.assertGreater(abs(pair_pressure), 0.)

        self.s.part[0].exclusions = [1]
        # Force and energy should not be changed by the exclusion
        self.s.integrator.run(0)
        self.assertAlmostEqual(self.s.part[0].f[0], pair_force, places=7)
        self.assertAlmostEqual(self.s.part[1].f[0], -pair_force, places=7)
        self.assertAlmostEqual(self.s.analysis.energy()[('coulomb', 0)], pair_energy, places=7)
        self.assertAlmostEqual(self.s.analysis.pressure()[('coulomb', 0)], pair_pressure, places=7)

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
