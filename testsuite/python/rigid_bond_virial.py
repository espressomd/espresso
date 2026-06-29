#
# Copyright (C) 2013-2026 The ESPResSo project
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

@utx.skipIfMissingFeatures("BOND_CONSTRAINT")
class RigidBondVirialTest(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4

    def setUp(self):
        self.system.time_step = 1e-2

    def tearDown(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def _make_dimer(self, m1=1.0, m2=1.0, v=1.0):
        x1 = 5. - m2 / (m1 + m2)
        x2 = 5. + m1 / (m1 + m2)
        v1 = (0., v, 0.)
        v2 = (0., -(m1 / m2) * v, 0.)
        bond = espressomd.interactions.RigidBond(r=1.0, ptol=1e-4, vtol=1e-4)
        self.system.bonded_inter.add(bond)
        p1 = self.system.part.add(pos=[x1, 5.0, 5.0], v=v1, mass=m1)
        p2 = self.system.part.add(pos=[x2, 5.0, 5.0], v=v2, mass=m2)
        p2.add_bond((bond, p1))
        return p1, p2

    #  Rotation test: after one step, constraint virial = F_centripetal·d
    #  F_centripetal = m v²/r = 1·1/0.5 = 2,  W = -virial/(3V) = -2/(3V)               
    def _virial_in_rotation(self, set_integrator):
        V = self.system.volume()
        set_integrator()
        self._make_dimer()
        self.system.integrator.run(1)
        v_p = np.trace(self.system.analysis.pressure_tensor()['bonded']) / 3.
        v_theory = -2.0 / (3. * V)
        self.assertAlmostEqual(v_p, v_theory, delta=0.01 * abs(v_theory))
        v_xx = self.system.analysis.pressure_tensor()['bonded'][0,0]
        v_yy = self.system.analysis.pressure_tensor()['bonded'][1,1]
        v_zz = self.system.analysis.pressure_tensor()['bonded'][2,2]
        self.assertAlmostEqual(v_xx, -2.0 / V, delta=0.01 * abs(2.0 / V))
        self.assertAlmostEqual(v_yy, 0.0, delta=1e-8)
        self.assertAlmostEqual(v_zz, 0.0, delta=1e-8)

    def test_virial_in_rotation_vv(self):
        """VV: constraint virial of a rotating rigid dimer matches centripetal theory."""
        self._virial_in_rotation(self.system.integrator.set_vv)

    def test_virial_in_rotation_se(self):
        """SE: constraint virial of a rotating rigid dimer matches centripetal theory."""
        self._virial_in_rotation(self.system.integrator.set_symplectic_euler)

    def _virial_unequal_masses(self, set_integrator):
        m1 = 1.0
        m2 = 2.0
        d = 1.0
        v = 1.0 # velocity
        V = self.system.volume()

        set_integrator()
        self._make_dimer(m1=m1, m2=m2, v=v)
        self.system.integrator.run(1)

        mu = m1 * m2 / (m1 + m2)               # reduced mass = 2/3
        omega = v * (m1 + m2) / (m2 * d)       # |v_rel| / d = 3/2
        v_theory = -mu * omega**2 * d**2 / (3. * V)   # = -1/(2V)

        pt = self.system.analysis.pressure_tensor()['bonded']
        v_p = np.trace(pt) / 3.
        self.assertAlmostEqual(v_p, v_theory, delta=0.01 * abs(v_theory))
        # Bond is along x: all virial goes into xx, none into yy or zz
        self.assertAlmostEqual(pt[0, 0], 3. * v_theory, delta=0.01 * abs(3. * v_theory))
        self.assertAlmostEqual(pt[1, 1], 0., delta=1e-8)
        self.assertAlmostEqual(pt[2, 2], 0., delta=1e-8)

    def test_virial_unequal_masses_vv(self):
        """VV: constraint virial with m1!=m2 matches centripetal theory."""
        self._virial_unequal_masses(self.system.integrator.set_vv)

    def test_virial_unequal_masses_se(self):
        """SE: constraint virial with m1!=m2 matches centripetal theory."""
        self._virial_unequal_masses(self.system.integrator.set_symplectic_euler)

    #  Langevin consistency: mean of (P_bond + P_kin) = kT/V,
    #  std matches analytic fluctuation formula.
    def _virial_consistency(self, set_integrator, noise_prefactor):
        kT = 1.0
        gamma = 1.0
        mass = 1.
        V = self.system.volume()
        dt = self.system.time_step
        self.system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        set_integrator()

        std_theory = ((6 * kT**2 + noise_prefactor * gamma * mass * kT / dt) / (9 * V**2))**0.5

        self._make_dimer(v=0.0)

        self.system.integrator.run(1000)   # equilibrate
        n_loop = 2000
        n_steps = 100
        virial = []
        for _ in range(n_loop):
            self.system.integrator.run(n_steps)
            v_p = np.trace(self.system.analysis.pressure_tensor()['bonded']) / 3.
            v_k = np.trace(self.system.analysis.pressure_tensor()['kinetic']) / 3.
            virial.append(v_p + v_k)

        rigid_p = np.mean(virial)
        rigid_std = np.std(virial)
        self.assertAlmostEqual(rigid_p, 1. / V, delta=2.*std_theory/n_loop**0.5)
        self.assertAlmostEqual(rigid_std, std_theory, delta=0.02*std_theory)

    def test_virial_consistency_vv(self):
        """VV+Langevin: rigid bond virial satisfies equipartition and fluctuation formula."""
        self._virial_consistency(self.system.integrator.set_vv, 1)

    def test_virial_consistency_se(self):
        """SE+Langevin: rigid bond virial satisfies equipartition and fluctuation formula."""
        self._virial_consistency(self.system.integrator.set_symplectic_euler, 4)


if __name__ == "__main__":
    ut.main()

