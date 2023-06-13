#
# Copyright (C) 2021-2023 The ESPResSo project
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
import espressomd.lb
import numpy as np


@utx.skipIfMissingFeatures("WALBERLA")
class Test(ut.TestCase):

    """
    Basic tests of the block forest.

    """
    system = espressomd.System(box_l=[12., 4., 4.])

    def test_interface(self):
        LatticeWalberla = espressomd.lb.LatticeWalberla

        # check getters
        for n_ghost_layers in range(10):
            obj = LatticeWalberla(agrid=1., n_ghost_layers=n_ghost_layers)
            self.assertEqual(obj.n_ghost_layers, n_ghost_layers)
        for agrid in (0.5, 1., 2.):
            obj = LatticeWalberla(agrid=agrid, n_ghost_layers=1)
            self.assertEqual(obj.agrid, agrid)
            target_shape = np.asarray(self.system.box_l, dtype=int) / obj.agrid
            np.testing.assert_array_equal(obj.shape, target_shape)

        # check exception mechanism
        obj = LatticeWalberla(agrid=1., n_ghost_layers=1)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'agrid' is read-only"):
            obj.agrid = 2.
        with self.assertRaisesRegex(RuntimeError, "Parameter 'n_ghost_layers' is read-only"):
            obj.n_ghost_layers = 2
        with self.assertRaisesRegex(RuntimeError, "Parameter 'n_ghost_layers' is missing"):
            LatticeWalberla(agrid=1.)
        with self.assertRaisesRegex(ValueError, "Parameter 'n_ghost_layers' must be >= 0"):
            LatticeWalberla(agrid=1., n_ghost_layers=-1)
        with self.assertRaisesRegex(ValueError, "Parameter 'agrid' must be > 0"):
            LatticeWalberla(agrid=0., n_ghost_layers=1)
        with self.assertRaisesRegex(ValueError, "Parameter 'agrid' must be > 0"):
            LatticeWalberla(agrid=-1., n_ghost_layers=1)
        with self.assertRaisesRegex(ValueError, "Parameter 'shape' must be derived from espressomd.shapes.Shape"):
            obj = LatticeWalberla(agrid=1., n_ghost_layers=1)
            next(obj.get_node_indices_inside_shape(10))


if __name__ == "__main__":
    ut.main()
