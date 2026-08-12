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

    # Rotation test: after one step, constraint virial = F_centripetal * d,
    # F_centripetal = m v^2/r = 1*1/0.5 = 2, W = -virial/(3V) = -2/(3V)
    def _virial_in_rotation(self, set_integrator):
        V = self.system.volume()
        set_integrator()
        self._make_dimer()
        self.system.integrator.run(1)
        v_p = np.trace(self.system.analysis.pressure_tensor()["bonded"]) / 3.
        v_theory = -2.0 / (3. * V)
        self.assertAlmostEqual(v_p, v_theory, delta=0.01 * abs(v_theory))
        v_xx = self.system.analysis.pressure_tensor()["bonded"][0, 0]
        v_yy = self.system.analysis.pressure_tensor()["bonded"][1, 1]
        v_zz = self.system.analysis.pressure_tensor()["bonded"][2, 2]
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
        d = 1.0  # distance between particle
        v = 1.0  # velocity
        V = self.system.volume()

        set_integrator()
        self._make_dimer(m1=m1, m2=m2, v=v)
        self.system.integrator.run(1)

        mu = m1 * m2 / (m1 + m2)               # reduced mass = 2/3
        omega = v * (m1 + m2) / (m2 * d)       # |v_rel| / d = 3/2
        v_theory = -mu * omega**2 * d**2 / (3. * V)   # = -1/(2V)

        pt = self.system.analysis.pressure_tensor()["bonded"]
        v_p = np.trace(pt) / 3.
        self.assertAlmostEqual(v_p, v_theory, delta=0.01 * abs(v_theory))
        # Bond is along x: all virial goes into xx, none into yy or zz
        self.assertAlmostEqual(pt[0, 0], 3. * v_theory,
                               delta=0.01 * abs(3. * v_theory))
        self.assertAlmostEqual(pt[1, 1], 0., delta=1e-8)
        self.assertAlmostEqual(pt[2, 2], 0., delta=1e-8)

    def test_virial_unequal_masses_vv(self):
        """VV: constraint virial with m1!=m2 matches centripetal theory."""
        self._virial_unequal_masses(self.system.integrator.set_vv)

    def test_virial_unequal_masses_se(self):
        """SE: constraint virial with m1!=m2 matches centripetal theory."""
        self._virial_unequal_masses(
            self.system.integrator.set_symplectic_euler)

    def _make_chain(self, masses, bond_lengths, omega=1.0):
        """
        Build a rigid, collinear 3-particle chain 0-1-2 (two RigidBonds,
        particle 1 shared by both) set up in rigid rotation about its
        center of mass, so the net constraint force on each particle is
        exactly the centripetal force F_i = -m_i * omega^2 * x_i.
        """
        l1, l2 = bond_lengths
        x1 = 0.0
        x0 = x1 - l1
        x2 = x1 + l2
        com = (masses[0] * x0 + masses[1] * x1 + masses[2] * x2) / sum(masses)
        x = [x0 - com, x1 - com, x2 - com]

        parts = []
        for m, xi in zip(masses, x):
            v = (0., omega * xi, 0.)
            p = self.system.part.add(pos=[5. + xi, 5.0, 5.0], v=v, mass=m)
            parts.append(p)

        bond01 = espressomd.interactions.RigidBond(r=l1, ptol=1e-9, vtol=1e-9)
        bond12 = espressomd.interactions.RigidBond(r=l2, ptol=1e-9, vtol=1e-9)
        self.system.bonded_inter.add(bond01)
        self.system.bonded_inter.add(bond12)
        parts[1].add_bond((bond01, parts[0]))
        parts[2].add_bond((bond12, parts[1]))
        return parts, x

    def _virial_chain(self, set_integrator, masses, bond_lengths, omega=1.0):
        V = self.system.volume()
        set_integrator()
        _, x = self._make_chain(masses, bond_lengths, omega=omega)
        self.system.integrator.run(1)

        # Centripetal theory, generalized from the dimer case: each particle's
        # net constraint force is F_i = -m_i*omega^2*x_i, so the constraint
        # virial is W_xx = sum_i x_i * F_i (bond is along x: all virial is
        # in xx, none in yy or zz).
        w_xx = sum(-m * omega**2 * xi**2 for m, xi in zip(masses, x))
        v_theory = w_xx / (3. * V)

        pt = self.system.analysis.pressure_tensor()["bonded"]
        v_p = np.trace(pt) / 3.
        self.assertAlmostEqual(v_p, v_theory, delta=0.01 * abs(v_theory))
        self.assertAlmostEqual(pt[0, 0], w_xx / V,
                               delta=0.01 * abs(w_xx / V))
        self.assertAlmostEqual(pt[1, 1], 0., delta=1e-8)
        self.assertAlmostEqual(pt[2, 2], 0., delta=1e-8)

    def test_virial_chain_symmetric_vv(self):
        """
        VV: constraint virial of a rotating rigid 3-particle chain
        (equal masses/bond lengths, shared middle bond) matches
        centripetal theory.
        """
        self._virial_chain(self.system.integrator.set_vv,
                           masses=[1.0, 1.0, 1.0], bond_lengths=(1.0, 1.0))

    def test_virial_chain_symmetric_se(self):
        """
        SE: constraint virial of a rotating rigid 3-particle chain
        (equal masses/bond lengths, shared middle bond) matches
        centripetal theory.
        """
        self._virial_chain(self.system.integrator.set_symplectic_euler,
                           masses=[1.0, 1.0, 1.0], bond_lengths=(1.0, 1.0))

    def test_virial_chain_unequal_masses_vv(self):
        """
        VV: constraint virial of a rigid 3-particle chain with unequal masses
        and bond lengths matches centripetal theory.
        """
        self._virial_chain(self.system.integrator.set_vv,
                           masses=[2.0, 1.0, 3.0], bond_lengths=(1.0, 1.5))

    def test_virial_chain_unequal_masses_se(self):
        """
        SE: constraint virial of a rigid 3-particle chain with unequal masses
        and bond lengths matches centripetal theory.
        """
        self._virial_chain(self.system.integrator.set_symplectic_euler,
                           masses=[2.0, 1.0, 3.0], bond_lengths=(1.0, 1.5))


if __name__ == "__main__":
    ut.main()
