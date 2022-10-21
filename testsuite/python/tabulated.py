#
# Copyright (C) 2013-2022 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.interactions
import numpy as np


class TabulatedTest(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.min_ = 1.
        self.max_ = 2.
        self.dx = (self.max_ - self.min_) / 99.
        self.force = 5 + np.arange(0, 100) * 2.3 * self.dx
        self.energy = 5 - np.arange(0, 100) * 2.3 * self.dx

        self.system.part.add(type=0, pos=[5., 5., 5.0])
        self.system.part.add(type=0, pos=[5., 5., 5.5])

    def tearDown(self):
        self.system.part.clear()

    def check(self):
        p0, p1 = self.system.part.all()
        # Below cutoff
        np.testing.assert_allclose(np.copy(self.system.part.all().f), 0.0)

        for z in np.linspace(0, self.max_ - self.min_, 200, endpoint=False):
            p1.pos = [5., 5., 6. + z]
            self.system.integrator.run(0)
            np.testing.assert_allclose(
                np.copy(p0.f), [0., 0., -(5. + z * 2.3)])
            np.testing.assert_allclose(np.copy(p0.f), -np.copy(p1.f))
            self.assertAlmostEqual(
                self.system.analysis.energy()['total'], 5. - z * 2.3)

    @utx.skipIfMissingFeatures("TABULATED")
    def test_non_bonded(self):
        self.system.non_bonded_inter[0, 0].tabulated.set_params(
            min=self.min_, max=self.max_, energy=self.energy, force=self.force)
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].tabulated.cutoff, self.max_)

        params = self.system.non_bonded_inter[0, 0].tabulated.get_params()
        np.testing.assert_allclose(params['force'], self.force)
        np.testing.assert_allclose(params['energy'], self.energy)
        self.assertAlmostEqual(params['min'], self.min_)
        self.assertAlmostEqual(params['max'], self.max_)

        self.check()

        self.system.non_bonded_inter[0, 0].tabulated.set_params(
            min=-1, max=-1, energy=[], force=[])
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].tabulated.cutoff, -1.)

        with self.assertRaisesRegex(ValueError, "TabulatedPotential parameter 'max' must be larger than or equal to parameter 'min'"):
            espressomd.interactions.TabulatedNonBonded(
                min=1., max=0., energy=[0.], force=[0.])
        with self.assertRaisesRegex(ValueError, "TabulatedPotential parameter 'force' must contain 1 element"):
            espressomd.interactions.TabulatedNonBonded(
                min=1., max=1., energy=[0., 0.], force=[0., 0.])
        with self.assertRaisesRegex(ValueError, "TabulatedPotential parameter 'force' must contain at least 1 element"):
            espressomd.interactions.TabulatedNonBonded(
                min=1., max=2., energy=[], force=[])
        with self.assertRaisesRegex(ValueError, "TabulatedPotential parameter 'force' must have the same size as parameter 'energy'"):
            espressomd.interactions.TabulatedNonBonded(
                min=1., max=2., energy=[0.], force=[0., 0.])

    @utx.skipIfMissingFeatures("TABULATED")
    def test_bonded(self):
        tb = espressomd.interactions.TabulatedDistance(
            min=self.min_, max=self.max_, energy=self.energy, force=self.force)
        self.system.bonded_inter.add(tb)

        np.testing.assert_allclose(tb.params['force'], self.force)
        np.testing.assert_allclose(tb.params['energy'], self.energy)
        self.assertAlmostEqual(tb.params['min'], self.min_)
        self.assertAlmostEqual(tb.params['max'], self.max_)

        p0, p1 = self.system.part.all()
        p0.add_bond((tb, p1))
        self.check()

        # make bond too short: potential becomes constant
        for z in np.linspace(0.1, 1., 9, endpoint=False):
            p1.pos = [5., 5., 5. + z]
            self.system.integrator.run(0)
            np.testing.assert_allclose(np.copy(p0.f), [0., 0., -5.])
            np.testing.assert_allclose(np.copy(p0.f), -np.copy(p1.f))
            self.assertAlmostEqual(self.system.analysis.energy()['total'], 5.)

        # break bond
        p1.pos = [5., 5., 6. + self.max_ + 0.1]
        with self.assertRaisesRegex(Exception, "bond broken between particles 0, 1"):
            self.system.analysis.energy()
        with self.assertRaisesRegex(Exception, "bond broken between particles 0, 1"):
            self.system.integrator.run(0)


if __name__ == "__main__":
    ut.main()
