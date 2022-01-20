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
import espressomd
from espressomd.bond_breakage import BreakageSpec
from espressomd.interactions import HarmonicBond


@utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
class BondBreakage(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        system = cls.system = espressomd.System(box_l=[10] * 3)
        system.cell_system.skin = 0.4
        system.time_step = 0.01
        system.min_global_cut = 2

        pos1 = system.box_l / 2 - 0.5
        pos2 = system.box_l / 2 + 0.5
        cls.p1 = system.part.add(pos=pos1)
        cls.p2 = system.part.add(pos=pos2)

        cls.p1v = system.part.add(pos=pos1)
        cls.p1v.vs_auto_relate_to(cls.p1)

        cls.p2v = system.part.add(pos=pos2)
        cls.p2v.vs_auto_relate_to(cls.p2)

        cls.h1 = HarmonicBond(k=1, r_0=0)
        cls.h2 = HarmonicBond(k=1, r_0=0)
        system.bonded_inter.add(cls.h1)
        system.bonded_inter.add(cls.h2)

    def tearDown(self):
        self.system.bond_breakage.clear()

    def test_00_interface(self):
        self.assertEqual(len(self.system.bond_breakage), 0)

        spec2 = BreakageSpec(breakage_length=1.2, action_type=1)
        spec4 = BreakageSpec(breakage_length=0.2, action_type=2)
        self.system.bond_breakage[2] = spec2
        self.system.bond_breakage[4] = spec4
        self.assertEqual(self.system.bond_breakage[2], spec2)
        self.assertEqual(self.system.bond_breakage[4], spec4)
        self.assertEqual(len(self.system.bond_breakage), 2)
        self.assertEqual(sorted(self.system.bond_breakage.keys()), [2, 4])
        self.assertEqual(
            sorted(self.system.bond_breakage.items()),
            [(2, spec2), (4, spec4)])

        self.system.bond_breakage.clear()
        self.assertEqual(len(self.system.bond_breakage), 0)
        self.assertEqual(self.system.bond_breakage.keys(), [])
        with self.assertRaisesRegex(ValueError, "Key has to be of type int"):
            self.system.bond_breakage[self.h1]
        with self.assertRaisesRegex(ValueError, "Key has to be of type int"):
            self.system.bond_breakage.remove(self.h1)
        with self.assertRaisesRegex(RuntimeError, "Inserting breakage spec without a bond type is not permitted"):
            self.system.bond_breakage.call_method("insert", object=spec2)

    def test_ignore(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1._bond_id] = BreakageSpec(
            breakage_length=2, action_type=1)

        self.p1.bonds = ((self.h1, self.p2))
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ((self.h1, self.p2.id),))

        self.p2.bonds = [(self.h1, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p2.bonds, ((self.h1, self.p1.id),))

        # Different bond type
        system.bond_breakage[self.h1._bond_id] = BreakageSpec(
            breakage_length=0.2, action_type=1)
        self.p1.bonds = [(self.h2, self.p2)]
        self.p2.bonds = [(self.h2, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ((self.h2, self.p2.id),))
        self.assertEqual(self.p2.bonds, ((self.h2, self.p1.id),))

    def test_delete_bond(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1._bond_id] = BreakageSpec(
            breakage_length=0, action_type=1)

        self.p1.bonds = [(self.h1, self.p2)]
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ())

        self.p2.bonds = [(self.h1, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p2.bonds, ())

    def test_revert_bind_at_point_of_collision(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1._bond_id] = BreakageSpec(
            breakage_length=0.5, action_type=2)

        self.p1.bonds = [(self.h2, self.p2)]
        self.p1v.bonds = [(self.h1, self.p2v)]
        system.integrator.run(1)
        self.assertEqual(self.p1v.bonds, ())
        self.assertEqual(self.p1.bonds, ())

        self.p2.bonds = [(self.h2, self.p1)]
        self.p1v.bonds = [(self.h1, self.p2v)]
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ())
        self.assertEqual(self.p1v.bonds, ())

    def test_exceptions(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h2._bond_id] = BreakageSpec(
            breakage_length=0.5, action_type=2)

        self.p1.bonds = [(self.h2, self.p2)]
        self.p1v.bonds = [(self.h1, self.p2v)]
        with self.assertRaisesRegex(Exception, "The REVERT_BIND_AT_POINT_OF_COLLISION bond breakage action has to be configured for the bond on the virtual site"):
            system.integrator.run(1)


if __name__ == "__main__":
    ut.main()
