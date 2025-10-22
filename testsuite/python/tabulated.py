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
import contextlib
with contextlib.suppress(ImportError):
    import sympy as sp


@utx.skipIfMissingFeatures("TABULATED")
class TabulatedTest(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.min_ = 1.
        self.max_ = 2.

        self.system.part.add(type=0, pos=[5., 5., 5.0])
        self.system.part.add(type=0, pos=[5., 5., 5.5])

    def tearDown(self):
        self.system.part.clear()
        self.system.non_bonded_inter.reset()
        self.system.bonded_inter.clear()

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

    def test_non_bonded(self):
        dx = (self.max_ - self.min_) / 99.
        ref_force = 5. + np.arange(0, 100) * 2.3 * dx
        ref_energy = 5. - np.arange(0, 100) * 2.3 * dx
        self.system.non_bonded_inter[0, 0].tabulated.set_params(
            min=self.min_, max=self.max_, energy=ref_energy, force=ref_force)
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].tabulated.cutoff, self.max_)

        params = self.system.non_bonded_inter[0, 0].tabulated.get_params()
        np.testing.assert_allclose(params['force'], ref_force)
        np.testing.assert_allclose(params['energy'], ref_energy)
        self.assertAlmostEqual(params['min'], self.min_)
        self.assertAlmostEqual(params['max'], self.max_)

        self.check()

        self.system.non_bonded_inter[0, 0].tabulated.set_params(
            min=-1, max=-1, energy=[], force=[])
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].tabulated.cutoff, -1.)

    def test_non_bonded_exceptions(self):
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

    def test_bonded(self):
        dx = (self.max_ - self.min_) / 99.
        ref_force = 5. + np.arange(0, 100) * 2.3 * dx
        ref_energy = 5. - np.arange(0, 100) * 2.3 * dx
        tb = espressomd.interactions.TabulatedDistance(
            min=self.min_, max=self.max_, energy=ref_energy, force=ref_force)
        self.system.bonded_inter.add(tb)

        np.testing.assert_allclose(tb.params['force'], ref_force)
        np.testing.assert_allclose(tb.params['energy'], ref_energy)
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

    @utx.skipIfMissingModules("sympy")
    def test_tabulated_sympy_string(self):
        eps = 5.
        sig = 1.2
        steps = 100
        r = np.linspace(self.min_, self.max_, steps)
        ref_energy = 4 * eps * ((sig / r)**12 - (sig / r)**6)
        ref_force = -4 * eps * (-12 / r * (sig / r)**12 + 6 / r * (sig / r)**6)
        ia_params = self.system.non_bonded_inter[0, 0].tabulated
        ia_params.set_analytical(
            min=self.min_, max=self.max_, steps=steps, sig=sig, eps=eps,
            energy_expr="4*eps*((sig/r)**12-(sig/r)**6)")
        np.testing.assert_allclose(np.copy(ia_params.energy), ref_energy)
        np.testing.assert_allclose(np.copy(ia_params.force), ref_force)
        np.testing.assert_allclose(ia_params.cutoff, self.max_)

    @utx.skipIfMissingModules("sympy")
    def test_tabulated_sympy_expr(self):
        steps = 20
        x = sp.Symbol("x")
        fn = sp.Expr(x**2)
        r = np.linspace(self.min_, self.max_, steps)
        ref_energy = r**2
        ref_force = -2 * r
        ia_params = self.system.non_bonded_inter[0, 0].tabulated
        ia_params.set_analytical(
            min=self.min_, max=self.max_, steps=steps, dist_param=x,
            energy_expr=fn)
        np.testing.assert_allclose(np.copy(ia_params.energy), ref_energy)
        np.testing.assert_allclose(np.copy(ia_params.force), ref_force)
        np.testing.assert_allclose(ia_params.cutoff, self.max_)

    @utx.skipIfMissingModules("sympy")
    def test_tabulated_sympy_expr_like(self):
        steps = 20
        eps = 5.
        x = sp.Symbol("x")
        e = sp.Symbol("eps")
        fn = e * x**2
        r = np.linspace(self.min_, self.max_, steps)
        ref_energy = eps * r**2
        ref_force = -2 * eps * r
        ref_fn_args = fn.args
        ia_params = self.system.non_bonded_inter[0, 0].tabulated
        ia_params.set_analytical(
            min=self.min_, max=self.max_, steps=steps, dist_param=x,
            energy_expr=fn, eps=5.)
        fn_args = fn.args
        np.testing.assert_allclose(np.copy(ia_params.energy), ref_energy)
        np.testing.assert_allclose(np.copy(ia_params.force), ref_force)
        np.testing.assert_allclose(ia_params.cutoff, self.max_)
        # check the input parameter wasn't modified
        self.assertEqual(fn_args, ref_fn_args)

    @utx.skipIfMissingModules("sympy")
    def test_tabulated_sympy_piecewise(self):
        steps = 20
        r = np.linspace(self.min_, self.max_, steps)
        # generate a curve that cuts off in the middle of the range
        midpoint = (self.min_ + self.max_) / 2.
        midpoint -= 1e-5 * midpoint  # offset to avoid floating-point ambiguity
        x = sp.Symbol("x")
        fn = sp.Piecewise((x**2, x < midpoint), (0., True))
        ref_energy = np.multiply(r**2, np.repeat([1, 0], steps // 2))
        ref_force = np.multiply(-2 * r, np.repeat([1, 0], steps // 2))
        ia_params = self.system.non_bonded_inter[0, 0].tabulated
        ia_params.set_analytical(
            min=self.min_, max=self.max_, steps=steps, dist_param=x,
            energy_expr=fn)
        np.testing.assert_allclose(np.copy(ia_params.energy), ref_energy)
        np.testing.assert_allclose(np.copy(ia_params.force), ref_force)
        np.testing.assert_allclose(ia_params.cutoff, self.max_)

    @utx.skipIfMissingModules("sympy")
    def test_tabulated_sympy_exceptions(self):
        ia_params = self.system.non_bonded_inter[0, 0].tabulated
        with self.assertRaisesRegex(TypeError, "Parameter 'dist_param' isn't sympy-compatible"):
            ia_params.set_analytical(
                min=0., max=1., steps=2, dist_param=None, energy_expr="x**2")
        with self.assertRaisesRegex(TypeError, "Parameter 'energy_expr' isn't sympy-compatible"):
            ia_params.set_analytical(
                min=0., max=1., steps=2, eps=1., energy_expr=None)
        with self.assertRaisesRegex(ValueError, "Parameter 'energy_expr' isn't a function of 'r'"):
            ia_params.set_analytical(
                min=0., max=1., steps=2, sig=1., energy_expr="(sig/theta)**12")
        with self.assertRaisesRegex(ValueError, "Parameter 'eps' isn't a constant in 'energy_expr'"):
            ia_params.set_analytical(
                min=0., max=1., steps=2, eps=1., energy_expr="(sig/r)**12")


if __name__ == "__main__":
    ut.main()
