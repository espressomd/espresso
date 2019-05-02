# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
from espressomd import lb
from espressomd import utils


class ArrayLockedTest(ut.TestCase):

    def test_locked_operators(self):
        v = utils.array_locked([1., 2., 3.])

        with self.assertRaises(ValueError):
            v[0] = 0
        with self.assertRaises(ValueError):
            v += [1, 1, 1]
        with self.assertRaises(ValueError):
            v -= [1, 1, 1]
        with self.assertRaises(ValueError):
            v *= [1, 1, 1]
        with self.assertRaises(ValueError):
            v /= [1, 1, 1]
        with self.assertRaises(ValueError):
            v //= [1, 1, 1]
        with self.assertRaises(ValueError):
            v %= [1, 1, 1]
        with self.assertRaises(ValueError):
            v **= [1, 1, 1]
        with self.assertRaises(ValueError):
            v <<= [1, 1, 1]
        with self.assertRaises(ValueError):
            v >>= [1, 1, 1]
        with self.assertRaises(ValueError):
            v &= [1, 1, 1]
        with self.assertRaises(ValueError):
            v |= [1, 1, 1]
        with self.assertRaises(ValueError):
            v ^= [1, 1, 1]

    def test_unlocked_operators(self):
        v = utils.array_locked([1, 2, 3])
        w = utils.array_locked([4, 5, 6])
        add = v + w
        sub = v - w

        self.assertTrue(isinstance(add, np.ndarray))
        self.assertTrue(isinstance(sub, np.ndarray))

        self.assertTrue(add.flags.writeable)
        self.assertTrue(sub.flags.writeable)

        np.testing.assert_array_equal(add, np.add(np.copy(v), np.copy(w)))
        np.testing.assert_array_equal(sub, np.subtract(np.copy(v), np.copy(w)))
        np.testing.assert_array_equal(sub, -(w - v))

    def test_copy_is_writeable(self):
        v = np.copy(utils.array_locked([1, 2, 3]))
        self.assertTrue(v.flags.writeable)

    def test_setter(self):
        v = utils.array_locked([1, 2, 3])
        v = [4, 5, 6]
        np.testing.assert_array_equal(v, [4, 5, 6])


class ArrayPropertyTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.box_l = [12.0, 12.0, 12.0]
    system.time_step = 0.01
    system.cell_system.skin = 0.01
    system.part.add(pos=[0, 0, 0])
    lbf = lb.LBFluid(agrid=0.5, dens=1, visc=1, tau=0.01)
    system.actors.add(lbf)

    def locked_operators(self, v):
        with self.assertRaises(ValueError):
            v[0] = 0
        with self.assertRaises(ValueError):
            v += [1, 1, 1]
        with self.assertRaises(ValueError):
            v -= [1, 1, 1]
        with self.assertRaises(ValueError):
            v *= [1, 1, 1]
        with self.assertRaises(ValueError):
            v /= [1, 1, 1]
        with self.assertRaises(ValueError):
            v //= [1, 1, 1]
        with self.assertRaises(ValueError):
            v %= [1, 1, 1]
        with self.assertRaises(ValueError):
            v **= [1, 1, 1]
        with self.assertRaises(ValueError):
            v <<= [1, 1, 1]
        with self.assertRaises(ValueError):
            v >>= [1, 1, 1]
        with self.assertRaises(ValueError):
            v &= [1, 1, 1]
        with self.assertRaises(ValueError):
            v |= [1, 1, 1]
        with self.assertRaises(ValueError):
            v ^= [1, 1, 1]

    def set_copy(self, v):
        cpy = np.copy(v)
        self.assertTrue(cpy.flags.writeable)

    def test_common(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].pos)
        self.locked_operators(self.system.part[0].v)
        self.locked_operators(self.system.part[0].f)
        self.locked_operators(self.system.part[0].pos_folded)

        # System
        self.locked_operators(self.system.box_l)

        # Check (allowed) setter
        # Particle
        self.system.part[0].pos = [2, 2, 2]
        self.assertTrue((self.system.part[0].pos == [2, 2, 2]).all())

        self.system.part[0].v = [2, 2, 2]
        self.assertTrue((self.system.part[0].v == [2, 2, 2]).all())

        self.system.part[0].f = [2, 2, 2]
        self.assertTrue((self.system.part[0].f == [2, 2, 2]).all())

        # System
        self.system.box_l = [2, 2, 2]
        self.assertTrue((self.system.box_l == [2, 2, 2]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].pos)
        self.set_copy(self.system.part[0].pos)
        self.set_copy(self.system.part[0].v)
        self.set_copy(self.system.part[0].f)
        self.set_copy(self.system.part[0].pos_folded)

        # System
        self.set_copy(self.system.box_l)

    @ut.skipIf(not espressomd.has_features(["ROTATION"]),
               "Features not available, skipping test!")
    def test_rotation(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].omega_lab)
        self.locked_operators(self.system.part[0].quat)
        self.locked_operators(self.system.part[0].rotation)
        self.locked_operators(self.system.part[0].omega_body)
        self.locked_operators(self.system.part[0].torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.locked_operators(self.system.part[0].ext_torque)

        # Check (allowed) setter
        # Particle
        self.system.part[0].quat = [0.5, 0.5, 0.5, 0.5]
        self.assertTrue(
            (self.system.part[0].quat == [0.5, 0.5, 0.5, 0.5]).all())

        self.system.part[0].omega_lab = [2, 2, 2]
        self.assertTrue((self.system.part[0].omega_lab == [2, 2, 2]).all())

        self.system.part[0].rotation = [1, 1, 1]
        self.assertTrue((self.system.part[0].rotation == [1, 1, 1]).all())

        self.system.part[0].omega_body = [2, 2, 2]
        self.assertTrue((self.system.part[0].omega_body == [2, 2, 2]).all())

        self.system.part[0].torque_lab = [2, 2, 2]
        self.assertTrue((self.system.part[0].torque_lab == [2, 2, 2]).all())

        if espressomd.has_features("EXTERNAL_FORCES"):
            self.system.part[0].ext_torque = [2, 2, 2]
            self.assertTrue(
                (self.system.part[0].ext_torque == [2, 2, 2]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].omega_lab)
        self.set_copy(self.system.part[0].quat)
        self.set_copy(self.system.part[0].rotation)
        self.set_copy(self.system.part[0].omega_body)
        self.set_copy(self.system.part[0].torque_lab)
        if espressomd.has_features("EXTERNAL_FORCES"):
            self.set_copy(self.system.part[0].ext_torque)

    @ut.skipIf(not espressomd.has_features(["ROTATIONAL_INERTIA"]),
               "Features not available, skipping test!")
    def test_rotational_inertia(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].rinertia)

        # Check (allowed) setter
        # Particle
        self.system.part[0].rinertia = [2, 2, 2]
        self.assertTrue((self.system.part[0].rinertia == [2, 2, 2]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].rinertia)

    @ut.skipIf(not espressomd.has_features(["EXTERNAL_FORCES"]),
               "Features not available, skipping test!")
    def test_external_forces(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].ext_force)
        self.locked_operators(self.system.part[0].fix)

        # Check (allowed) setter
        # Particle
        self.system.part[0].ext_force = [2, 2, 2]
        self.assertTrue((self.system.part[0].ext_force == [2, 2, 2]).all())

        self.system.part[0].fix = [1, 1, 1]
        self.assertTrue((self.system.part[0].fix == [1, 1, 1]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].ext_force)
        self.set_copy(self.system.part[0].fix)

    @ut.skipIf(
        not espressomd.has_features(["ROTATION", "PARTICLE_ANISOTROPY"]),
           "Features not available, skipping test!")
    def test_rot_aniso(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].gamma_rot)

        # Check (allowed) setter
        # Particle
        self.system.part[0].gamma_rot = [2, 2, 2]
        self.assertTrue((self.system.part[0].gamma_rot == [2, 2, 2]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].gamma_rot)

    def test_lb(self):
        # Check for exception for various operators
        # LB
        self.locked_operators(self.lbf[0, 0, 0].velocity)
        self.locked_operators(self.lbf[0, 0, 0].stress)
        self.locked_operators(self.lbf[0, 0, 0].stress_neq)
        self.locked_operators(self.lbf[0, 0, 0].population)

    @ut.skipIf(not espressomd.has_features(["LANGEVIN_PER_PARTICLE",
                                            "PARTICLE_ANISOTROPY"]),
               "Features not available, skipping test!")
    def test_langevinpp_aniso(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].gamma)

        # Check (allowed) setter
        # Particle
        self.system.part[0].gamma = [2, 2, 2]
        self.assertTrue((self.system.part[0].gamma == [2, 2, 2]).all())

        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].gamma)

    @ut.skipIf(not espressomd.has_features(["DIPOLES"]),
               "Features not available, skipping test!")
    def test_dipoles(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].dip)

        # Check (allowed) setter
        # Particle
        self.system.part[0].dip = [2, 2, 2]
        np.testing.assert_allclose(
            [2, 2, 2], np.copy(self.system.part[0].dip), atol=1E-15)
        
        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].dip)

    @ut.skipIf(not espressomd.has_features(["EXCLUSIONS"]),
               "Features not available, skipping test!")
    def test_exclusions(self):

        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].exclusions)

    @ut.skipIf(not espressomd.has_features(["PARTIAL_PERIODIC"]),
               "Features not available, skipping test!")
    def test_partial_periodic(self):

        # Check for exception for various operators
        # System
        self.locked_operators(self.system.periodicity)

        # Check (allowed) setter
        # System
        self.system.periodicity = [1, 0, 0]
        self.assertTrue((self.system.periodicity == [1, 0, 0]).all())

        # Check if copy is settable
        # System
        self.set_copy(self.system.periodicity)

if __name__ == "__main__":
    ut.main()
