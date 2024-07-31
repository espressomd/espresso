#
# Copyright (C) 2024 The ESPResSo project
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
import espressomd.plugins.ase
import numpy as np
import ase


class ASEInterfaceTest(ut.TestCase):

    system = espressomd.System(
        box_l=[10., 10., 10.], periodicity=[True, True, False])

    def setUp(self):
        self.system.part.add(pos=[0., 0., 0.], f=[1., -1., 0.], type=0)
        self.system.part.add(pos=[0., 0., 1.], f=[0., 12., 0.], type=1)
        self.system.ase = espressomd.plugins.ase.ASEInterface(
            type_mapping={0: "H", 1: "O"},
        )

    def tearDown(self):
        self.system.part.clear()

    def test_ase_get(self):
        """Test the ``ASEInterface.get()`` method."""
        # Create a simple ASE atoms object
        atoms = self.system.ase.get()
        self.assertIsInstance(atoms, ase.Atoms)
        self.assertEqual(set(atoms.get_chemical_symbols()), {"H", "O"})
        np.testing.assert_equal(atoms.pbc, np.copy(self.system.periodicity))
        np.testing.assert_allclose(atoms.cell, np.diag(self.system.box_l))
        np.testing.assert_allclose(atoms.get_positions(),
                                   [[0., 0., 0.], [0., 0., 1.]])
        np.testing.assert_allclose(atoms.get_forces(),
                                   [[1., -1., 0.], [0., 12., 0.]])

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_exceptions(self):
        p = self.system.part.add(pos=[0., 0., 0.], type=10)
        with self.assertRaisesRegex(RuntimeError, r"Particle types '\{10\}' haven't been registered"):
            self.system.ase.get()
        p.type = 1
        vs = self.system.part.add(pos=[0., 0., 0.], type=1)
        vs.vs_relative = (p.id, 0.01, (1., 0., 0., 0.))
        with self.assertRaisesRegex(RuntimeError, "ASE doesn't support virtual sites"):
            self.system.ase.get()


if __name__ == "__main__":
    ut.main()
