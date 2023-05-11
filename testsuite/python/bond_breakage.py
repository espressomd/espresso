#
# Copyright (C) 2013-2022 The ESPResSo project
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
import numpy as np
import espressomd.interactions
from espressomd.bond_breakage import BreakageSpec
from espressomd.interactions import HarmonicBond, AngleHarmonic


class BondBreakageCommon:
    system = espressomd.System(box_l=[10] * 3)
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    system.min_global_cut = 2


@utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
class BondBreakage(BondBreakageCommon, ut.TestCase):

    @classmethod
    def setUpClass(cls):
        pos1 = cls.system.box_l / 2 - 0.5
        pos2 = cls.system.box_l / 2 + 0.5
        pos3 = cls.system.box_l / 2 - 1.5
        cls.p1 = cls.system.part.add(pos=pos1)
        cls.p2 = cls.system.part.add(pos=pos2)
        cls.p3 = cls.system.part.add(pos=pos3)

        cls.p1v = cls.system.part.add(pos=pos1)
        cls.p1v.vs_auto_relate_to(cls.p1)

        cls.p2v = cls.system.part.add(pos=pos2)
        cls.p2v.vs_auto_relate_to(cls.p2)

        cls.h1 = HarmonicBond(k=1, r_0=0)
        cls.h2 = HarmonicBond(k=1, r_0=0)
        cls.angle = AngleHarmonic(bend=1, phi0=np.pi)
        cls.system.bonded_inter.add(cls.h1)
        cls.system.bonded_inter.add(cls.h2)
        cls.system.bonded_inter.add(cls.angle)

    @classmethod
    def tearDownClass(cls):
        cls.system.part.clear()
        cls.system.bonded_inter.clear()

    def tearDown(self):
        self.system.bond_breakage.clear()

    def test_00_interface(self):
        self.assertEqual(len(self.system.bond_breakage), 0)

        spec2 = BreakageSpec(breakage_length=1.2, action_type="delete_bond")
        spec4 = BreakageSpec(breakage_length=0.2,
                             action_type="revert_bind_at_point_of_collision")
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
        with self.assertRaisesRegex(ValueError, "Key has to be of type 'int'"):
            self.system.bond_breakage[None]
        with self.assertRaisesRegex(ValueError, "Key has to be of type 'int'"):
            self.system.bond_breakage[None] = spec2
        with self.assertRaisesRegex(ValueError, "Key has to be of type 'int', got type 'double'"):
            self.system.bond_breakage.remove(1.)
        with self.assertRaisesRegex(ValueError, "Bond needs to be added to the system first"):
            self.system.bond_breakage[HarmonicBond(k=1, r_0=0)]
        with self.assertRaisesRegex(RuntimeError, "Inserting breakage spec without a bond type is not permitted"):
            self.system.bond_breakage.call_method("insert", object=spec2)

    def test_ignore(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1] = BreakageSpec(
            breakage_length=2, action_type="delete_bond")

        self.p1.bonds = ((self.h1, self.p2))
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ((self.h1, self.p2.id),))

        self.p2.bonds = [(self.h1, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p2.bonds, ((self.h1, self.p1.id),))

        # Different bond type
        system.bond_breakage[self.h1] = BreakageSpec(
            breakage_length=0.2, action_type="delete_bond")
        self.p1.bonds = [(self.h2, self.p2)]
        self.p2.bonds = [(self.h2, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ((self.h2, self.p2.id),))
        self.assertEqual(self.p2.bonds, ((self.h2, self.p1.id),))

    def test_delete_bond(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1] = BreakageSpec(
            breakage_length=0, action_type="delete_bond")

        self.p1.bonds = [(self.h1, self.p2)]
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ())

        self.p2.bonds = [(self.h1, self.p1)]
        system.integrator.run(1)
        self.assertEqual(self.p2.bonds, ())

    def test_delete_angle_bond(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.angle] = BreakageSpec(
            breakage_length=5, action_type="delete_bond")
        self.p1.bonds = ((self.angle, self.p2, self.p3),)
        bonds = self.p1.bonds
        system.integrator.run(1)
        # should still be there. Not bfeyond breakage disst
        self.assertEqual(self.p1.bonds, bonds)
        self.system.bond_breakage.clear()
        system.bond_breakage[self.angle] = BreakageSpec(
            breakage_length=1, action_type="delete_bond")
        system.integrator.run(1)
        self.assertEqual(self.p1.bonds, ())

    def test_revert_bind_at_point_of_collision_pair(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h1] = BreakageSpec(
            breakage_length=0.5, action_type="revert_bind_at_point_of_collision")

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

    def test_revert_bind_at_point_of_collision_angle(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.angle] = BreakageSpec(
            breakage_length=0.5, action_type="revert_bind_at_point_of_collision")

        self.p1.bonds = [(self.h2, self.p2)]
        self.p2.bonds = [(self.h2, self.p1)]
        self.p1v.bonds = [(self.angle, self.p1, self.p2)]
        self.p2v.bonds = [(self.angle, self.p1, self.p2)]
        system.integrator.run(1)
        self.assertEqual(self.p1v.bonds, ())
        self.assertEqual(self.p2v.bonds, ())
        self.assertEqual(self.p1.bonds, ())
        self.assertEqual(self.p2.bonds, ())

    def test_exceptions(self):
        system = self.system

        # Particles closer than cutoff
        system.bond_breakage[self.h2] = BreakageSpec(
            breakage_length=0.5, action_type="revert_bind_at_point_of_collision")

        self.p1.bonds = [(self.h2, self.p2)]
        self.p1v.bonds = [(self.h1, self.p2v)]
        with self.assertRaisesRegex(Exception, "The REVERT_BIND_AT_POINT_OF_COLLISION bond breakage action has to be configured for the bond on the virtual site"):
            system.integrator.run(1)

        self.system.bond_breakage.clear()

        # Particles closer than cutoff
        system.bond_breakage[self.angle] = BreakageSpec(
            breakage_length=0.5, action_type="revert_bind_at_point_of_collision")

        self.p1.bonds = [(self.h2, self.p2)]
        self.p2.bonds = [(self.h2, self.p1)]
        self.p1.bonds = [(self.angle, self.p1v, self.p2)]
        with self.assertRaisesRegex(Exception, "The REVERT_BIND_AT_POINT_OF_COLLISION bond breakage action has to be configured for the bond on the virtual site"):
            system.integrator.run(1)


@utx.skipIfMissingFeatures(["LENNARD_JONES", "COLLISION_DETECTION"])
class NetworkBreakage(BondBreakageCommon, ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.system.box_l = 3 * [20]
        cls.system.min_global_cut = 0.6
        cls.system.time_step = 0.01
        cls.system.cell_system.skin = 0.4

    def count_bonds(self, pairs):
        bonds_count = 0
        for pair in pairs:
            for bond in self.system.part.by_id(pair[0]).bonds:
                if bond[1] == pair[1]:
                    bonds_count += 1
            for bond in self.system.part.by_id(pair[1]).bonds:
                if bond[1] == pair[0]:
                    bonds_count += 1
        return bonds_count

    def setUp(self):

        box_vol = self.system.box_l[0]**3.
        phi = 0.4

        r = 1.
        solid_vol = phi * box_vol
        part_vol = 4 / 3 * np.pi * r**3
        part_num = int(solid_vol / part_vol)

        np.random.seed(seed=678)
        for i in range(part_num):
            pos = np.random.rand(3) * self.system.box_l[0]
            self.system.part.add(pos=pos)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=1., epsilon=1., cutoff=2**(1 / 6), shift='auto')
        self.system.integrator.set_steepest_descent(f_max=0,
                                                    gamma=0.1,
                                                    max_displacement=0.1)
        self.system.integrator.run(100)
        self.system.integrator.set_vv()

        for i in range(part_num):
            self.system.part.by_id(i).fix = [True, True, True]

        self.system.thermostat.set_langevin(kT=0.0, gamma=1.0, seed=41)

    def tearDown(self):
        self.system.part.clear()
        self.system.bonded_inter.clear()
        self.system.thermostat.turn_off()

    @utx.skipIfMissingFeatures(["COLLISION_DETECTION"])
    def test_center_bonds(self):

        harm = espressomd.interactions.HarmonicBond(k=1.0, r_0=0.0, r_cut=5)
        self.system.bonded_inter.add(harm)

        crit = 2**(1 / 6) * 2.

        self.system.collision_detection.set_params(mode="bind_centers",
                                                   distance=2**(1 / 6) * 2.2, bond_centers=harm)
        self.system.integrator.run(1)

        self.system.collision_detection.set_params(mode="off")
        self.system.bond_breakage[harm] = BreakageSpec(
            breakage_length=crit, action_type="delete_bond")
        self.system.integrator.run(1)

        bonds_dist = 0
        pairs = self.system.cell_system.get_pairs(crit, types=[0])
        for pair in pairs:
            dist = self.system.distance(
                self.system.part.by_id(pair[0]),
                self.system.part.by_id(pair[1]))
            if dist <= crit:
                bonds_dist += 1

        bonds_count = self.count_bonds(pairs)
        np.testing.assert_equal(bonds_dist, bonds_count)

    @utx.skipIfMissingFeatures(
        ["VIRTUAL_SITES_RELATIVE", "COLLISION_DETECTION"])
    def test_vs_bonds(self):

        harm = espressomd.interactions.HarmonicBond(k=1.0, r_0=0.0, r_cut=5)
        virt = espressomd.interactions.Virtual()
        self.system.bonded_inter.add(harm)
        self.system.bonded_inter.add(virt)

        crit = 2**(1 / 6) * 1.5
        crit_vs = 2**(1 / 6) * 1 / 3 * 1.2

        self.system.collision_detection.set_params(mode="bind_at_point_of_collision",
                                                   distance=crit, bond_centers=virt, bond_vs=harm,
                                                   part_type_vs=1, vs_placement=1 / 3)
        self.system.integrator.run(1)

        self.system.collision_detection.set_params(mode="off")
        self.system.bond_breakage[harm] = BreakageSpec(
            breakage_length=crit_vs, action_type="revert_bind_at_point_of_collision")
        self.system.integrator.run(1)

        bonds_dist = 0
        pairs = self.system.cell_system.get_pairs(
            2**(1 / 6) * 2 / 3, types=[1])

        for pair in pairs:
            r1 = self.system.part.by_id(pair[0]).vs_relative[0]
            r2 = self.system.part.by_id(pair[1]).vs_relative[0]
            dist = self.system.distance(
                self.system.part.by_id(r1),
                self.system.part.by_id(r2))
            dist_vs = self.system.distance(
                self.system.part.by_id(pair[0]),
                self.system.part.by_id(pair[1]))
            if dist_vs <= crit_vs:
                if dist <= crit:
                    if dist > 0.0:
                        bonds_dist += 1

        bonds_count = self.count_bonds(pairs)

        np.testing.assert_equal(bonds_dist, bonds_count)


if __name__ == "__main__":
    ut.main()
