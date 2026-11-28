#
# Copyright (C) 2026 The ESPResSo project
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


@utx.skipIfMissingFeatures(["BOND_CONSTRAINT", "MASS"])
class RigidBondVirialStatsTest(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4

    def setUp(self):
        self.system.time_step = 1e-2

    def tearDown(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def _make_dimer(self):
        bond = espressomd.interactions.RigidBond(r=1.0, ptol=1e-4, vtol=1e-4)
        self.system.bonded_inter.add(bond)
        p1 = self.system.part.add(pos=[4.5, 5.0, 5.0], v=[0., 0., 0.], mass=1.)
        p2 = self.system.part.add(pos=[5.5, 5.0, 5.0], v=[0., 0., 0.], mass=1.)
        p2.add_bond((bond, p1))
        return p1, p2

    # Langevin consistency: mean of (P_bond + P_kin) = kT/V,
    # standard deviation matches analytic fluctuation formula.
    def _virial_consistency(self, set_integrator, noise_prefactor):
        kT = 1.0
        gamma = 1.0
        mass = 1.
        V = self.system.volume()
        dt = self.system.time_step
        self.system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
        set_integrator()

        std_theory = ((6 * kT**2 + noise_prefactor * gamma *
                      mass * kT / dt) / (9 * V**2))**0.5

        self._make_dimer()

        self.system.integrator.run(1000)  # equilibrate
        n_loop = 2000
        n_steps = 100
        virial = []
        for _ in range(n_loop):
            self.system.integrator.run(n_steps)
            v_p = np.trace(
                self.system.analysis.pressure_tensor()["bonded"]) / 3.
            v_k = np.trace(self.system.analysis.pressure_tensor()[
                           "kinetic"]) / 3.
            virial.append(v_p + v_k)

        rigid_p = np.mean(virial)
        rigid_std = np.std(virial)
        self.assertAlmostEqual(
            rigid_p, 1. / V, delta=2. * std_theory / n_loop**0.5)
        self.assertAlmostEqual(rigid_std, std_theory, delta=0.02 * std_theory)

    def test_virial_consistency_vv(self):
        """VV+Langevin: rigid bond virial satisfies equipartition and fluctuation formula."""
        self._virial_consistency(self.system.integrator.set_vv, 1)

    def test_virial_consistency_se(self):
        """SE+Langevin: rigid bond virial satisfies equipartition and fluctuation formula."""
        self._virial_consistency(
            self.system.integrator.set_symplectic_euler, 4)


if __name__ == "__main__":
    ut.main()
