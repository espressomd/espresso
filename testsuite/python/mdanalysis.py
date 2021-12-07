#
# Copyright (C) 2019-2020 The ESPResSo project
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

"""
Testmodule for the MDAnalysis interface.
"""
import espressomd
import espressomd.interactions
import numpy as np
import unittest as ut
import unittest_decorators as utx
try:
    import MDAnalysis as mda
    import espressomd.MDA_ESP
    skipIfMissingPythonPackage = utx.no_skip
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module MDAnalysis not available, skipping test!")


@skipIfMissingPythonPackage
class TestMDAnalysis(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.001
    system.cell_system.skin = 0.1

    for i in range(10):
        system.part.add(id=i, pos=[i, i % 2, 0], v=[0, i, -i], f=[1, 2 * i, 0],
                        type=i % 2, q=i % 3 - 1)

    bond = espressomd.interactions.HarmonicBond(k=1., r_0=1.2, r_cut=2.0)
    angle = espressomd.interactions.AngleCosine(bend=1., phi0=2 * np.pi / 3)
    dihe = espressomd.interactions.Dihedral(bend=1., mult=2, phase=np.pi / 3)
    system.bonded_inter.add(bond)
    system.bonded_inter.add(angle)
    system.bonded_inter.add(dihe)
    system.part.by_id(1).add_bond((bond, 0))
    system.part.by_id(3).add_bond((angle, 2, 4))
    system.part.by_id(6).add_bond((dihe, 5, 7, 8))

    def test_universe(self):
        system = self.system
        partcls = system.part.all()
        eos = espressomd.MDA_ESP.Stream(system)
        u = mda.Universe(eos.topology, eos.trajectory)
        # check atoms
        self.assertEqual(len(u.atoms), 10)
        np.testing.assert_equal(u.atoms.ids, np.arange(10) + 1)
        np.testing.assert_equal(u.atoms.types, 5 * ['T0', 'T1'])
        np.testing.assert_almost_equal(
            u.atoms.charges, partcls.q, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.positions, partcls.pos, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.velocities, partcls.v, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.forces, partcls.f, decimal=6)
        # check bonds
        self.assertEqual(len(u.bonds), 1)
        self.assertEqual(len(u.angles), 1)
        self.assertEqual(len(u.dihedrals), 1)
        np.testing.assert_equal(u.bonds[0].atoms.ids, [1, 2])
        np.testing.assert_equal(u.angles[0].atoms.ids, [3, 4, 5])
        np.testing.assert_equal(u.dihedrals[0].atoms.ids, [6, 7, 8, 9])
        self.assertEqual(len(u.atoms[9].bonds), 0)
        self.assertEqual(len(u.atoms[9].angles), 0)
        self.assertEqual(len(u.atoms[9].dihedrals), 0)


if __name__ == "__main__":
    ut.main()
