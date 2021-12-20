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
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.electrostatics


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class ElectrostaticInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System(box_l=[20., 20., 20.])

    def setUp(self):
        self.system.time_step = 0.01

        self.p0 = self.system.part.add(pos=(9.0, 2.0, 2.0), q=1)
        self.p1 = self.system.part.add(pos=(11.0, 2.0, 2.0), q=-1)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def calc_dh_potential(self, r, dh_params):
        kT = 1.0
        q1, q2 = self.system.part.all().q
        u = np.zeros_like(r)
        # r<r_cut
        i = np.where(r < dh_params['r_cut'])[0]
        u[i] = dh_params['prefactor'] * kT * q1 * \
            q2 * np.exp(-dh_params['kappa'] * r[i]) / r[i]
        return u

    def calc_rf_potential(self, r, rf_params):
        """Calculates the potential of the ReactionField coulomb method"""

        kT = 1.0

        q1, q2 = self.system.part.all().q
        epsilon1 = rf_params['epsilon1']
        epsilon2 = rf_params['epsilon2']
        kappa = rf_params['kappa']
        r_cut = rf_params['r_cut']
        # prefactor calculation
        B = (2 * (epsilon1 - epsilon2) * (1 + kappa * r_cut) -
             epsilon2 * kappa * kappa * r_cut * r_cut) / \
            ((epsilon1 + 2 * epsilon2) * (1 + kappa * r_cut) +
             epsilon2 * kappa * kappa * r_cut * r_cut)
        offset = (1. - B / 2.) / r_cut
        u = np.zeros_like(r)

        # r<r_cut
        i = np.where(r < rf_params['r_cut'])[0]
        u[i] = rf_params['prefactor'] * kT * q1 * q2 * \
            ((1. / r[i] - B * np.square(r[i]) / (2. * r_cut**3)) - offset)
        return u

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m(self):
        prefactor = 1.1
        box_vol = self.system.volume()
        dip = np.copy(self.p0.q * self.p0.pos + self.p1.q * self.p1.pos)
        p3m_params = {'accuracy': 1e-7,
                      'mesh': [22, 22, 22],
                      'cao': 7,
                      'r_cut': 8.906249999999998,
                      'alpha': 0.387611049779351}

        # reference values for energy and force calculated for prefactor = 1
        ref_energy = -0.501062398379 * prefactor
        ref_force1 = [0.248921612 * prefactor, 0, 0]
        ref_force2 = [-ref_force1[0], 0, 0]

        # check metallic case
        p3m = espressomd.electrostatics.P3M(
            prefactor=prefactor, epsilon='metallic', tune=False, **p3m_params)
        self.system.actors.add(p3m)
        self.system.integrator.run(0, recalc_forces=True)
        p3m_energy = self.system.analysis.energy()['coulomb']
        tol = 1e-5
        np.testing.assert_allclose(p3m_energy, ref_energy, atol=tol)
        np.testing.assert_allclose(np.copy(self.p0.f), ref_force1, atol=tol)
        np.testing.assert_allclose(np.copy(self.p1.f), ref_force2, atol=tol)

        # keep current values as reference to check for P3M dipole correction
        ref_energy_metallic = self.system.analysis.energy()['coulomb']
        ref_forces_metallic = np.copy(self.system.part.all().f)
        self.system.actors.remove(p3m)

        # check non-metallic case
        tol = 1e-10
        for epsilon in np.power(10., np.arange(-4, 5)):
            dipole_correction = 4 * np.pi / box_vol / (1 + 2 * epsilon)
            energy_correction = dipole_correction * np.linalg.norm(dip)**2
            forces_correction = np.outer(
                [self.p0.q, self.p1.q], dipole_correction * dip)
            ref_energy = ref_energy_metallic + prefactor * energy_correction
            ref_forces = ref_forces_metallic - prefactor * forces_correction
            p3m = espressomd.electrostatics.P3M(
                prefactor=prefactor, epsilon=epsilon, tune=False, **p3m_params)
            self.system.actors.add(p3m)
            self.system.integrator.run(0, recalc_forces=True)
            p3m_forces = np.array([self.p0.f, self.p1.f])
            p3m_energy = self.system.analysis.energy()['coulomb']
            np.testing.assert_allclose(p3m_energy, ref_energy, atol=tol)
            np.testing.assert_allclose(p3m_forces, ref_forces, atol=tol)
            self.system.actors.remove(p3m)

    def test_dh(self):
        dh_params = dict(prefactor=1.2, kappa=0.8, r_cut=2.0)
        dh = espressomd.electrostatics.DH(
            prefactor=dh_params['prefactor'],
            kappa=dh_params['kappa'],
            r_cut=dh_params['r_cut'])

        self.system.actors.add(dh)
        dr = 0.001
        r = np.arange(.5, 1.01 * dh_params['r_cut'], dr)
        u_dh = self.calc_dh_potential(r, dh_params)
        f_dh = -np.gradient(u_dh, dr)
        # zero the discontinuity, and re-evaluate the derivative as a backwards
        # difference
        i_cut = np.argmin((dh_params['r_cut'] - r)**2)
        f_dh[i_cut] = 0
        f_dh[i_cut - 1] = (u_dh[i_cut - 2] - u_dh[i_cut - 1]) / dr

        u_dh_core = np.zeros_like(r)
        f_dh_core = np.zeros_like(r)

        for i, ri in enumerate(r):
            self.p1.pos = self.p0.pos + [ri, 0, 0]
            self.system.integrator.run(0)
            u_dh_core[i] = self.system.analysis.energy()['coulomb']
            f_dh_core[i] = self.p0.f[0]

        np.testing.assert_allclose(u_dh_core, u_dh, atol=1e-7)
        np.testing.assert_allclose(f_dh_core, -f_dh, atol=1e-2)

    def test_dh_pure_coulomb(self):
        dh_params = dict(prefactor=1.2, kappa=0.0, r_cut=2.0)
        dh = espressomd.electrostatics.DH(
            prefactor=dh_params['prefactor'],
            kappa=dh_params['kappa'],
            r_cut=dh_params['r_cut'])

        self.system.actors.add(dh)
        dr = 0.001
        r = np.arange(.5, 1.01 * dh_params['r_cut'], dr)
        u_dh = self.calc_dh_potential(r, dh_params)
        f_dh = u_dh / r

        u_dh_core = np.zeros_like(r)
        f_dh_core = np.zeros_like(r)

        for i, ri in enumerate(r):
            self.p1.pos = self.p0.pos + [ri, 0, 0]
            self.system.integrator.run(0)
            u_dh_core[i] = self.system.analysis.energy()['coulomb']
            f_dh_core[i] = self.p0.f[0]

        np.testing.assert_allclose(u_dh_core, u_dh, atol=1e-7)
        np.testing.assert_allclose(f_dh_core, -f_dh, atol=1e-7)

    def test_dh_exceptions(self):
        dh = espressomd.electrostatics.DH(prefactor=-1.0, kappa=1.0, r_cut=1.0)
        with self.assertRaisesRegex(ValueError, 'Coulomb prefactor has to be >= 0'):
            self.system.actors.add(dh)
        self.system.actors.clear()
        dh = espressomd.electrostatics.DH(prefactor=1.0, kappa=-1.0, r_cut=1.0)
        with self.assertRaisesRegex(ValueError, 'kappa should be a non-negative number'):
            self.system.actors.add(dh)
        self.system.actors.clear()
        dh = espressomd.electrostatics.DH(prefactor=1.0, kappa=1.0, r_cut=-1.0)
        with self.assertRaisesRegex(ValueError, 'r_cut should be a non-negative number'):
            self.system.actors.add(dh)

    def test_rf(self):
        """Tests the ReactionField coulomb interaction by comparing the
           potential and force against the analytic values"""

        rf_params = dict(prefactor=1.0,
                         kappa=2.0,
                         epsilon1=1.0,
                         epsilon2=2.0,
                         r_cut=2.0)
        rf = espressomd.electrostatics.ReactionField(
            prefactor=rf_params['prefactor'],
            kappa=rf_params['kappa'],
            epsilon1=rf_params['epsilon1'],
            epsilon2=rf_params['epsilon2'],
            r_cut=rf_params['r_cut'])
        self.system.actors.add(rf)

        dr = 0.001
        r = np.arange(.5, 1.01 * rf_params['r_cut'], dr)

        u_rf = self.calc_rf_potential(r, rf_params)
        f_rf = -np.gradient(u_rf, dr)

        # zero the discontinuity, and re-evaluate the derivative as a backwards
        # difference
        i_cut = np.argmin((rf_params['r_cut'] - r)**2)
        f_rf[i_cut] = 0
        f_rf[i_cut - 1] = (u_rf[i_cut - 2] - u_rf[i_cut - 1]) / dr

        u_rf_core = np.zeros_like(r)
        f_rf_core = np.zeros_like(r)

        for i, ri in enumerate(r):
            self.p1.pos = self.p0.pos + [ri, 0, 0]
            self.system.integrator.run(0)
            u_rf_core[i] = self.system.analysis.energy()['coulomb']
            f_rf_core[i] = self.p0.f[0]

        np.testing.assert_allclose(u_rf_core, u_rf, atol=1e-7)
        np.testing.assert_allclose(f_rf_core, -f_rf, atol=1e-2)

    def test_rf_exceptions(self):
        params = dict(kappa=1.0, epsilon1=1.0, epsilon2=2.0, r_cut=1.0)
        for key in params:
            invalid_params = {**params, 'prefactor': 1.0, key: -1.0}
            rf = espressomd.electrostatics.ReactionField(**invalid_params)
            with self.assertRaisesRegex(ValueError, f'{key} should be a non-negative number'):
                self.system.actors.add(rf)
            self.system.actors.clear()


if __name__ == "__main__":
    ut.main()
