#
# Copyright (C) 2021 The ESPResSo project
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


@utx.skipIfMissingFeatures("WALBERLA")
class TestLattice(ut.TestCase):

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

        # check exception mechanism
        obj = LatticeWalberla(agrid=1., n_ghost_layers=1)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'agrid' is read-only"):
            obj.agrid = 2.
        with self.assertRaisesRegex(RuntimeError, "Parameter 'n_ghost_layers' is read-only"):
            obj.n_ghost_layers = 2
        with self.assertRaisesRegex(ValueError, 'The following keys have to be given as keyword arguments'):
            LatticeWalberla(agrid=1.)
        with self.assertRaises(ValueError):
            LatticeWalberla(agrid=1., n_ghost_layers=-1)
        with self.assertRaises(ValueError):
            LatticeWalberla(agrid=0., n_ghost_layers=1)
        with self.assertRaises(ValueError):
            LatticeWalberla(agrid=-1., n_ghost_layers=1)


if __name__ == "__main__":
    ut.main()
