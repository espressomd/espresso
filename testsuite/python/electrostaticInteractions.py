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
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [20, 20, 20]
        self.system.time_step = 0.01

        if not self.system.part.exists(0):
            self.system.part.add(id=0, pos=(1.0, 2.0, 2.0), q=1)
        if not self.system.part.exists(1):
            self.system.part.add(
                id=1, pos=(3.0, 2.0, 2.0), q=-1)

    def calc_dh_potential(self, r, df_params):
        kT = 1.0
        q1 = self.system.part[0].q
        q2 = self.system.part[1].q
        u = np.zeros_like(r)
        # r<r_cut
        i = np.where(r < df_params['r_cut'])[0]
        u[i] = df_params['prefactor'] * kT * q1 * \
            q2 * np.exp(-df_params['kappa'] * r[i]) / r[i]
        return u

    def calc_rf_potential(self, r, rf_params):
        """Calculates the potential of the ReactionField coulomb method"""

        kT = 1.0

        q1 = self.system.part[0].q
        q2 = self.system.part[1].q
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
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos = [3.0, 2.0, 2.0]
        # results, reference values for energy and force only calculated for
        # prefactor = 1
        p3m_energy = -0.501062398379 * prefactor
        p3m_force = 2.48921612e-01 * prefactor
        p3m = espressomd.electrostatics.P3M(prefactor=prefactor,
                                            accuracy=9.910945054074526e-08,
                                            mesh=[22, 22, 22],
                                            cao=7,
                                            r_cut=8.906249999999998,
                                            alpha=0.387611049779351,
                                            tune=False)
        self.system.actors.add(p3m)
        self.assertAlmostEqual(self.system.analysis.energy()['coulomb'],
                               p3m_energy, places=5)
        # need to update forces
        self.system.integrator.run(0)
        np.testing.assert_allclose(np.copy(self.system.part[0].f),
                                   [p3m_force, 0, 0], atol=1E-4)
        np.testing.assert_allclose(np.copy(self.system.part[1].f),
                                   [-p3m_force, 0, 0], atol=1E-5)
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
        # need to update forces
        for i, ri in enumerate(r):
            self.system.part[1].pos = self.system.part[0].pos + [ri, 0, 0]
            self.system.integrator.run(0)
            u_dh_core[i] = self.system.analysis.energy()['coulomb']
            f_dh_core[i] = self.system.part[0].f[0]

        np.testing.assert_allclose(u_dh_core, u_dh, atol=1e-7)
        np.testing.assert_allclose(f_dh_core, -f_dh, atol=1e-2)
        self.system.actors.remove(dh)

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
            self.system.part[1].pos = self.system.part[0].pos + [ri, 0, 0]
            self.system.integrator.run(0)
            u_rf_core[i] = self.system.analysis.energy()['coulomb']
            f_rf_core[i] = self.system.part[0].f[0]

        np.testing.assert_allclose(u_rf_core, u_rf, atol=1e-7)
        np.testing.assert_allclose(f_rf_core, -f_rf, atol=1e-2)
        self.system.actors.remove(rf)


if __name__ == "__main__":
    ut.main()
