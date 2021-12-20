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
# Integration test for exclusions
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(['EXCLUSIONS'])
class Exclusions(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.part.clear()
        self.system.box_l = 3 * [10]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01

    def test_add_remove(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0])
        self.system.part.add(id=1, pos=[0, 0, 0])
        self.system.part.add(id=2, pos=[0, 0, 0])

        p0.add_exclusion(1)
        p0.add_exclusion(2)
        self.assertEqual(list(p0.exclusions), [1, 2])
        p0.delete_exclusion(1)
        self.assertEqual(list(p0.exclusions), [2])
        p0.delete_exclusion(2)
        self.assertEqual(list(p0.exclusions), [])

    def test_transfer(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0], v=[1., 1., 1])
        self.system.part.add(id=1, pos=[0, 0, 0])
        self.system.part.add(id=2, pos=[0, 0, 0])
        self.system.part.add(id=3, pos=[0, 0, 0])

        p0.exclusions = [1, 2, 3]

        for _ in range(15):
            self.system.integrator.run(100)
            self.assertEqual(list(p0.exclusions), [1, 2, 3])

    @utx.skipIfMissingFeatures(['LENNARD_JONES'])
    def test_particle_property(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=2., cutoff=1.5, shift=0.0)

        p0 = self.system.part.add(id=0, pos=[0, 0, 0], type=0)
        p1 = self.system.part.add(id=1, pos=[1, 0, 0], type=0)

        pair_energy = self.system.analysis.energy()['total']
        self.assertGreater(pair_energy, 0.)

        pair_pressure = self.system.analysis.pressure()['total']
        self.assertGreater(pair_pressure, 0.)

        self.system.integrator.run(0)
        pair_force = p0.f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)

        p2 = self.system.part.add(id=2, pos=[2, 0, 0], type=0)
        self.system.integrator.run(0)
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               2 * pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               2 * pair_pressure)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = [0, 2]
        self.system.integrator.run(0)
        self.assertAlmostEqual(self.system.analysis.energy()['total'], 0)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'], 0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], 0, places=7)
        self.assertAlmostEqual(p2.f[0], 0, places=7)

        p1.exclusions = [0]
        self.assertAlmostEqual(
            self.system.analysis.energy()['total'],
            pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], pair_force, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = []
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               2 * pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               2 * pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], pair_force, places=7)
        self.assertAlmostEqual(p1.f[0], 0, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = [0]
        self.assertAlmostEqual(
            self.system.analysis.energy()['total'],
            pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], pair_force, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

    @utx.skipIfMissingFeatures(['P3M'])
    def test_electrostatics_not_excluded(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0], type=0, q=+1.)
        p1 = self.system.part.add(id=1, pos=[1, 0, 0], type=0, q=-1.)

        # Small alpha means large short-range contribution
        p3m = espressomd.electrostatics.P3M(
            prefactor=1, r_cut=3.0, accuracy=1e-3, mesh=32, cao=7, alpha=0.1,
            tune=False)
        self.system.actors.add(p3m)

        # Only short-range part of the coulomb energy
        pair_energy = self.system.analysis.energy()[('coulomb', 0)]
        self.assertGreater(abs(pair_energy), 0.)

        self.system.integrator.run(0)
        pair_force = p0.f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)

        pair_pressure = self.system.analysis.pressure()[('coulomb', 0)]
        self.assertGreater(abs(pair_pressure), 0.)

        p0.exclusions = [1]
        # Force and energy should not be changed by the exclusion
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], pair_force, places=7)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)
        self.assertAlmostEqual(self.system.analysis.energy()[('coulomb', 0)],
                               pair_energy, places=7)
        self.assertAlmostEqual(self.system.analysis.pressure()[('coulomb', 0)],
                               pair_pressure, places=7)


if __name__ == "__main__":
    ut.main()
