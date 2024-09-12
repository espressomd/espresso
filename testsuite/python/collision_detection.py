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
import espressomd.interactions
import espressomd.propagation
import numpy as np
import random


@utx.skipIfMissingFeatures("COLLISION_DETECTION")
class CollisionDetection(ut.TestCase):

    """Tests functionality of the collision detection / dynamic binding"""

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(seed=42)

    bond_center = espressomd.interactions.HarmonicBond(k=5000, r_0=0.1)
    bond_vs = espressomd.interactions.HarmonicBond(k=25000, r_0=0.02)
    bond_pair = espressomd.interactions.HarmonicBond(k=100, r_0=0.1)
    bond_angle_vs = espressomd.interactions.AngleHarmonic(
        bend=0., phi0=np.pi / 3.)
    system.bonded_inter.add(bond_center)
    system.bonded_inter.add(bond_vs)
    system.bonded_inter.add(bond_pair)
    system.bonded_inter.add(bond_angle_vs)
    time_step = 0.001
    system.time_step = time_step
    system.cell_system.skin = 0.05
    system.min_global_cut = 0.112

    part_type_to_attach_vs_to = 0
    part_type_vs = 1
    part_type_to_be_glued = 2
    part_type_after_glueing = 3
    other_type = 5

    def tearDown(self):
        self.system.time_step = self.time_step

    def get_state_set_state_consistency(self):
        state = self.system.collision_detection.get_params()
        self.system.collision_detection.set_params(**state)
        self.assertEqual(state, self.system.collision_detection.get_params())

    def test_bind_centers(self):
        system = self.system
        # Check that it leaves particles alone, when off
        system.collision_detection.protocol = espressomd.collision_detection.Off()

        system.part.clear()
        p0 = system.part.add(pos=(0, 0, 0), id=0)
        p1 = system.part.add(pos=(0.1, 0, 0), id=1)
        p2 = system.part.add(pos=(0.1, 0.3, 0), id=2)
        system.integrator.run(1)
        self.assertEqual(p0.bonds, ())
        self.assertEqual(p1.bonds, ())
        self.assertEqual(p2.bonds, ())
        # Check that it cannot be activated
        system.collision_detection.protocol = espressomd.collision_detection.BindCenters(
            distance=0.11, bond_centers=self.bond_center)
        self.get_state_set_state_consistency()
        system.integrator.run(1, recalc_forces=True)
        bond0 = ((system.bonded_inter[0], 1),)
        bond1 = ((system.bonded_inter[0], 0),)
        self.assertTrue(
            p0.bonds == bond0 or p1.bonds == bond1)
        self.assertEqual(p2.bonds, ())

        # Check that no additional bonds appear
        system.integrator.run(1)
        self.assertTrue(
            p0.bonds == bond0 or p1.bonds == bond1)
        self.assertEqual(p2.bonds, ())

        # Check turning it off
        system.collision_detection.protocol = espressomd.collision_detection.Off()
        self.get_state_set_state_consistency()
        self.assertIsInstance(
            self.system.collision_detection.protocol, espressomd.collision_detection.Off)

    def run_test_bind_at_point_of_collision_for_pos(self, *positions):
        system = self.system
        positions = list(positions)
        random.shuffle(positions)
        system.part.clear()
        # Place particle which should not take part in collisions
        p = system.part.add(pos=(0.1, 0.3, 0.))
        for pos in positions:
            p1 = system.part.add(pos=pos + (0., 0., 0.))
            p2 = system.part.add(pos=pos + (0.1, 0., 0.))
            assert system.distance(p1, p) >= 0.12 and system.distance(
                p2, p) >= 0.12, "Test particles too close to particle, which should not take part in collision"

        # 2 non-virtual + 2 virtual + one that doesn't take part
        expected_np = 4 * len(positions) + 1

        system.collision_detection.protocol = espressomd.collision_detection.BindAtPointOfCollision(
            bond_centers=self.bond_center, bond_vs=self.bond_vs, part_type_vs=1, vs_placement=0.4, distance=0.11)
        self.get_state_set_state_consistency()
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Integrate again and check that nothing has changed
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Check that nothing explodes when the particles are moved.
        # In particular for parallel simulations
        system.thermostat.set_langevin(kT=0, gamma=0.01, seed=42)
        system.part.all().v = [0.05, 0.01, 0.15]
        system.integrator.run(3000)
        self.verify_state_after_bind_at_poc(expected_np)

    def verify_state_after_bind_at_poc(self, expected_np):
        if self.system.collision_detection.protocol.bond_vs == self.bond_angle_vs:
            self.verify_state_after_bind_at_poc_triplet(expected_np)
        else:
            self.verify_state_after_bind_at_poc_pair(expected_np)

    def verify_state_after_bind_at_poc_pair(self, expected_np):
        system = self.system
        self.assertEqual(len(system.part), expected_np)
        Propagation = espressomd.propagation.Propagation

        # At the end of test, this list should be empty
        parts_not_accounted_for = list(range(expected_np))

        # We traverse particles. We look for a vs with a bond to find the other vs.
        # From the two vs we find the two non-virtual particles
        for p in system.part:
            # Skip non-virtual
            if not p.propagation & Propagation.TRANS_VS_RELATIVE:
                continue
            # Skip vs that doesn't have a bond
            if p.bonds == ():
                continue
            # Parse the bond
            self.assertEqual(len(p.bonds), 1)
            # Bond type
            self.assertEqual(p.bonds[0][0], self.bond_vs)
            # get partner
            p2 = system.part.by_id(p.bonds[0][1])
            # Is that really a vs
            self.assertTrue(p2.is_virtual())
            # Get base particles
            base_p1 = system.part.by_id(p.vs_relative[0])
            base_p2 = system.part.by_id(p2.vs_relative[0])
            # Take note of accounted-for particles
            for _p in (p, p2, base_p1, base_p2):
                parts_not_accounted_for.remove(_p.id)
            self.verify_bind_at_poc(base_p1, base_p2, p, p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for), 1)
        p = system.part.by_id(parts_not_accounted_for[0])
        self.assertFalse(p.propagation & Propagation.TRANS_VS_RELATIVE)
        self.assertEqual(p.bonds, ())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for, [])

    def verify_state_after_bind_at_poc_triplet(self, expected_np):
        system = self.system
        self.assertEqual(len(system.part), expected_np)

        # At the end of test, this list should be empty
        parts_not_accounted_for = list(range(expected_np))

        # We traverse particles. We look for a vs with a bond to find the other vs.
        # From the two vs we find the two non-virtual particles
        for p in system.part:
            # Skip non-virtual
            if not p.is_virtual():
                continue
            # Skip vs that doesn't have a bond
            if p.bonds == ():
                continue
            if p.id == 1:
                # Parse the bond
                self.assertEqual(len(p.bonds), 2)
                # Bond type
                self.assertEqual(p.bonds[0][0], self.bond_pair)
                self.assertEqual(p.bonds[1][0], self.bond_vs)
            else:
                # Parse the bond
                self.assertEqual(len(p.bonds), 1)
                # Bond type
                self.assertEqual(p.bonds[0][0], self.bond_angle_vs)
                # Is that really a vs
                self.assertTrue(p.is_virtual())
                # Get base particles
                base_p1 = system.part.by_id(p.bonds[0][1])
                base_p2 = system.part.by_id(p.bonds[0][2])
                # Take note of accounted-for particles
                for _p in (p, base_p1, base_p2):
                    if _p.id in parts_not_accounted_for:
                        parts_not_accounted_for.remove(_p.id)
        self.verify_bind_at_poc(system.part.by_id(1), system.part.by_id(2),
                                system.part.by_id(3), system.part.by_id(4))
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for), 1)
        p = system.part.by_id(parts_not_accounted_for[0])
        self.assertFalse(p.is_virtual())
        self.assertEqual(p.bonds, ())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for, [])

    def verify_bind_at_poc(self, p1, p2, vs1, vs2):
        system = self.system
        # Check for presence of vs
        # Check for bond between vs
        if self.system.collision_detection.protocol.bond_vs == self.bond_angle_vs:
            bond_p1 = ((self.bond_pair, p2.id), (self.bond_center, p2.id),)
            bond_p2 = ((self.bond_pair, p1.id), (self.bond_center, p1.id),)
            self.assertTrue(p1.bonds == bond_p1 or p2.bonds == bond_p2)
            bond_vs1 = ((self.bond_angle_vs, p1.id, p2.id),)
            bond_vs2 = ((self.bond_angle_vs, p2.id, p1.id),)
            self.assertTrue(vs1.bonds == bond_vs1 or vs1.bonds == bond_vs2)
            self.assertTrue(vs2.bonds == bond_vs1 or vs2.bonds == bond_vs2)
        else:
            bond_p1 = ((self.bond_center, p2.id),)
            bond_p2 = ((self.bond_center, p1.id),)
            self.assertTrue(p1.bonds == bond_p1 or p2.bonds == bond_p2)
            bond_vs1 = ((self.bond_vs, vs2.id),)
            bond_vs2 = ((self.bond_vs, vs1.id),)
            self.assertTrue(vs1.bonds == bond_vs1 or vs2.bonds == bond_vs2)

        # Vs properties
        self.assertTrue(vs1.is_virtual())
        self.assertTrue(vs2.is_virtual())

        # vs_relative properties
        seen = []
        for p in vs1, vs2:
            r = p.vs_relative
            rel_to = r[0]
            dist = r[1]
            # Vs is related to one of the particles
            self.assertIn(rel_to, (p1.id, p2.id))
            # The two vs relate to two different particles
            self.assertNotIn(rel_to, seen)
            seen.append(rel_to)

            # Check placement
            if rel_to == p1.id:
                dist_centers = np.copy(p2.pos - p1.pos)
            else:
                dist_centers = p1.pos - p2.pos
            expected_pos = system.part.by_id(rel_to).pos_folded + \
                system.collision_detection.protocol.vs_placement * dist_centers
            dist = expected_pos - p.pos_folded
            dist -= np.round(dist / system.box_l) * system.box_l
            self.assertLess(np.linalg.norm(dist), 1E-12)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_bind_at_point_of_collision(self):
        # Single collision head node
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0, 0, 0)))
        # Single collision, mixed
        self.run_test_bind_at_point_of_collision_for_pos(
            np.array((0.45, 0, 0)))
        # Single collision, non-head-node
        self.run_test_bind_at_point_of_collision_for_pos(np.array((0.7, 0, 0)))

        # head-node + mixed
        self.run_test_bind_at_point_of_collision_for_pos(
            np.array((0, 0, 0)), np.array((0.45, 0, 0)))
        # Mixed + other node
        self.run_test_bind_at_point_of_collision_for_pos(
            np.array((0.45, 0, 0)), np.array((0.7, 0, 0)))
        # Head + other
        self.run_test_bind_at_point_of_collision_for_pos(
            np.array((0.0, 0, 0)), np.array((0.7, 0, 0)))
        # Head + mixed + other
        self.run_test_bind_at_point_of_collision_for_pos(
            np.array((0.2, 0, 0)), np.array((0.95, 0, 0)), np.array((0.7, 0, 0)))

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_bind_at_point_of_collision_triplet(self):
        system = self.system
        positions = [np.array((0, 0, 0))]
        random.shuffle(positions)
        system.part.clear()
        # Place particle which should not take part in collisions
        p = system.part.add(pos=(0.1, 0.3, 0.))
        for pos in positions:
            p1 = system.part.add(pos=pos + (0., 0., 0.))
            p2 = system.part.add(pos=pos + (0.1, 0., 0.))
            p1.add_bond((self.bond_pair, p2))
            assert system.distance(p1, p) >= 0.12 and system.distance(
                p2, p) >= 0.12, "Test particles too close to particle, which should not take part in collision"

        # 2 non-virtual + 2 virtual + one that doesn't take part
        expected_np = 4 * len(positions) + 1

        system.collision_detection.protocol = espressomd.collision_detection.BindAtPointOfCollision(
            bond_centers=self.bond_center,
            bond_vs=self.bond_angle_vs, part_type_vs=1, vs_placement=0.4, distance=0.11)
        self.get_state_set_state_consistency()
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Integrate again and check that nothing has changed
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Check that nothing explodes when the particles are moved.
        # In particular for parallel simulations
        system.thermostat.set_langevin(kT=0, gamma=0.01, seed=42)
        system.part.all().v = [0.05, 0.01, 0.15]
        system.integrator.run(3000)
        self.verify_state_after_bind_at_poc(expected_np)

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "VIRTUAL_SITES_RELATIVE"])
    def test_bind_at_point_of_collision_random(self):
        """Integrate lj liquid and check that no double bonds are formed
           and the number of bonds fits the number of virtual sites

        """
        system = self.system
        system.part.clear()

        # Add randomly placed particles
        system.part.add(pos=np.random.random((200, 3)))

        # Setup Lennard-Jones
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=0.1, cutoff=2**(1. / 6) * 0.1, shift="auto")

        # Remove overlap between particles
        system.thermostat.turn_off()
        system.integrator.set_steepest_descent(
            f_max=0,
            gamma=1,
            max_displacement=0.001)
        while system.analysis.energy()["total"] > len(system.part):
            system.integrator.run(10)

        # Collision detection
        system.collision_detection.protocol = espressomd.collision_detection.BindAtPointOfCollision(
            bond_centers=self.bond_center, bond_vs=self.bond_vs, part_type_vs=1, vs_placement=0.4, distance=0.11)
        self.get_state_set_state_consistency()

        # Integrate lj liquid
        system.integrator.set_vv()
        system.integrator.run(5000)

        # Analysis
        virtual_sites = system.part.select(lambda p: p.is_virtual() == True)
        non_virtual = system.part.select(lambda p: p.is_virtual() == False)

        # Check bonds on non-virtual particles
        bonds = []
        for p in non_virtual:
            for bond in p.bonds:
                # Sort bond partners to make them unique independently of
                # which particle got the bond
                bonds.append(tuple(sorted([p.id, bond[1]])))

        # No duplicate bonds?
        self.assertEqual(len(bonds), len(set(bonds)))

        # 2 virtual sites per bond?
        self.assertEqual(2 * len(bonds), len(virtual_sites))

        # Find pairs of bonded virtual sites
        vs_pairs = []
        for p in virtual_sites:
            # 0 or 1 bond on vs?
            self.assertIn(len(p.bonds), [0, 1])

            if len(p.bonds) == 1:
                vs_pairs.append((p.id, p.bonds[0][1]))

        # Number of vs pairs = number of bonds?
        self.assertEqual(len(vs_pairs), len(bonds))

        # Check that vs pairs and bonds agree
        for vs_pair in vs_pairs:
            # Get corresponding non-virtual particles
            base_particles = tuple(sorted(
                [system.part.by_id(vs_pair[0]).vs_relative[0],
                 system.part.by_id(vs_pair[1]).vs_relative[0]]))

            # Is there a corresponding bond?
            self.assertIn(base_particles, bonds)

        # Tidy
        system.non_bonded_inter[0, 0].lennard_jones.deactivate()

    def run_test_glue_to_surface_for_pos(self, *positions):
        system = self.system
        positions = list(positions)
        random.shuffle(positions)
        system.part.clear()
        # Place particle which should not take part in collisions
        # In this case, it is skipped, because it is of the wrong type,
        # even if it is within range for a collision
        system.part.add(pos=positions[0], type=self.other_type)
        for pos in positions:
            # Since this is non-symmetric, we randomize order
            if np.random.random() > .5:
                system.part.add(
                    pos=pos + (0, 0, 0), type=self.part_type_to_attach_vs_to)
                system.part.add(
                    pos=pos + (0.1, 0, 0), type=self.part_type_to_be_glued)
            else:
                system.part.add(
                    pos=pos + (0.1, 0, 0), type=self.part_type_to_be_glued)
                system.part.add(
                    pos=pos + (0, 0, 0), type=self.part_type_to_attach_vs_to)

        # 2 non-virtual + 1 virtual + one that doesn't take part
        expected_np = 3 * len(positions) + 1

        system.collision_detection.protocol = espressomd.collision_detection.GlueToSurface(
            distance=0.11,
            distance_glued_particle_to_vs=0.02, bond_centers=self.bond_center,
            bond_vs=self.bond_vs, part_type_vs=self.part_type_vs,
            part_type_to_attach_vs_to=self.part_type_to_attach_vs_to,
            part_type_to_be_glued=self.part_type_to_be_glued,
            part_type_after_glueing=self.part_type_after_glueing)
        self.get_state_set_state_consistency()
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)

        # Integrate again and check that nothing has changed
        system.integrator.run(1, recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)

        # Check that nothing explodes, when the particles are moved.
        # In particular for parallel simulations
        system.thermostat.set_langevin(kT=0, gamma=0.01, seed=42)
        system.part.all().v = [0.05, 0.01, 0.15]
        system.integrator.run(3000)
        self.verify_state_after_glue_to_surface(expected_np)

    def verify_state_after_glue_to_surface(self, expected_np):
        system = self.system
        self.assertEqual(len(system.part), expected_np)
        Propagation = espressomd.propagation.Propagation

        # At the end of test, this list should be empty
        parts_not_accounted_for = list(range(expected_np))

        # We traverse particles. We look for a vs, get base particle from there
        # and partner particle via bonds
        for p in system.part:
            # Skip non-virtual
            if not p.propagation & Propagation.TRANS_VS_RELATIVE:
                continue
            # The vs shouldn't have bonds
            self.assertEqual(p.bonds, ())

            # Get base particles
            base_p = system.part.by_id(p.vs_relative[0])

            # Get bound particle
            # There is a bond between the base particle and the bound particle
            # but we have no guarantee, on where it is stored
            # 1. On the base particle of the vs
            p2 = None
            if len(base_p.bonds) == 1:
                self.assertEqual(base_p.bonds[0][0], self.bond_center)
                p2 = system.part.by_id(base_p.bonds[0][1])
            else:
                # We need to go through all particles to find it
                for candidate in system.part:
                    if candidate.id not in parts_not_accounted_for:
                        continue
                    if len(candidate.bonds) >= 1:
                        for b in candidate.bonds:
                            if b[0] == self.bond_center and b[1] == base_p.id:
                                p2 = candidate
                assert p2 is not None, "Bound particle not found"
            # Take note of accounted-for particles
            parts_not_accounted_for.remove(base_p.id)
            parts_not_accounted_for.remove(p.id)
            parts_not_accounted_for.remove(p2.id)
            self.verify_glue_to_surface_pair(base_p, p, p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for), 1)
        p = system.part.by_id(parts_not_accounted_for[0])
        self.assertFalse(p.propagation & Propagation.TRANS_VS_RELATIVE)
        self.assertEqual(p.type, self.other_type)
        self.assertEqual(p.bonds, ())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for, [])

    def verify_glue_to_surface_pair(self, base_p, vs, bound_p):
        system = self.system
        Propagation = espressomd.propagation.Propagation
        # Check all types
        self.assertEqual(base_p.type, self.part_type_to_attach_vs_to)
        self.assertEqual(vs.type, self.part_type_vs)
        self.assertEqual(bound_p.type, self.part_type_after_glueing)

        # Bound particle should have a bond to vs. It can additionally have a bond
        # to the base particle
        bond_to_vs_found = 0
        for b in bound_p.bonds:
            if b[0] == self.bond_vs:
                # bond to vs
                self.assertEqual(b, (self.bond_vs, vs.id))
                bond_to_vs_found += 1
        self.assertEqual(bond_to_vs_found, 1)
        # Vs should not have a bond
        self.assertEqual(vs.bonds, ())

        # Vs properties
        self.assertTrue(vs.propagation & Propagation.TRANS_VS_RELATIVE)
        self.assertEqual(vs.vs_relative[0], base_p.id)

        # Distance vs,bound_p
        self.assertAlmostEqual(system.distance(vs, bound_p), 0.02, places=3)
        self.assertAlmostEqual(system.distance(base_p, bound_p), 0.1, places=3)
        self.assertAlmostEqual(system.distance(base_p, vs), 0.08, places=3)

        # base_p,vs,bound_p on a line
        self.assertGreater(np.dot(system.distance_vec(base_p, vs), system.distance_vec(base_p, bound_p))
                           / system.distance(base_p, vs) / system.distance(base_p, bound_p), 0.99)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_glue_to_surface(self):
        # Single collision head node
        self.run_test_glue_to_surface_for_pos(np.array((0, 0, 0)))
        # Single collision, mixed
        self.run_test_glue_to_surface_for_pos(np.array((0.45, 0, 0)))
        # Single collision, non-head-node
        self.run_test_glue_to_surface_for_pos(np.array((0.7, 0, 0)))

        # head-node + mixed
        self.run_test_glue_to_surface_for_pos(
            np.array((0, 0, 0)), np.array((0.45, 0, 0)))
        # Mixed + other node
        self.run_test_glue_to_surface_for_pos(
            np.array((0.45, 0, 0)), np.array((0.7, 0, 0)))
        # Head + other
        self.run_test_glue_to_surface_for_pos(
            np.array((0.0, 0, 0)), np.array((0.7, 0, 0)))
        # Head + mixed + other
        self.run_test_glue_to_surface_for_pos(
            np.array((0.2, 0, 0)), np.array((0.95, 0, 0)), np.array((0.7, 0, 0)))

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "VIRTUAL_SITES_RELATIVE"])
    def test_glue_to_surface_random(self):
        """Integrate lj liquid and check that no double bonds are formed
           and the number of bonds fits the number of virtual sites

        """
        system = self.system
        system.part.clear()
        Propagation = espressomd.propagation.Propagation

        # Add randomly placed particles
        system.part.add(pos=np.random.random((300, 3)),
                        type=np.repeat([self.part_type_to_attach_vs_to,
                                        self.part_type_to_be_glued,
                                        self.other_type], [100, 100, 100]))

        # Setup Lennard-Jones
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=0.1, cutoff=2**(1. / 6) * 0.1, shift="auto")

        # Remove overlap between particles
        system.thermostat.turn_off()
        system.integrator.set_steepest_descent(
            f_max=0,
            gamma=1,
            max_displacement=0.001)
        while system.analysis.energy()["total"] > len(system.part):
            system.integrator.run(10)

        # Collision detection
        system.collision_detection.protocol = espressomd.collision_detection.GlueToSurface(
            distance=0.11,
            distance_glued_particle_to_vs=0.02, bond_centers=self.bond_center,
            bond_vs=self.bond_vs, part_type_vs=self.part_type_vs,
            part_type_to_attach_vs_to=self.part_type_to_attach_vs_to,
            part_type_to_be_glued=self.part_type_to_be_glued,
            part_type_after_glueing=self.part_type_after_glueing)
        self.get_state_set_state_consistency()

        # Integrate lj liquid
        system.integrator.set_vv()
        system.integrator.run(500)

        # Analysis
        virtual_sites = system.part.select(
            lambda p: p.propagation & Propagation.TRANS_VS_RELATIVE)
        non_virtual = system.part.select(lambda p: not (
            p.propagation & Propagation.TRANS_VS_RELATIVE))
        after_glueing = system.part.select(type=self.part_type_after_glueing)

        # One virtual site per glued particle?
        self.assertEqual(len(after_glueing), len(virtual_sites))

        # Check bonds on non-virtual particles
        bonds_centers = []
        bonds_virtual = []
        for p in non_virtual:
            # Inert particles should not have bonds
            if p.type == self.other_type:
                self.assertEqual(len(p.bonds), 0)

            # Particles that have not yet collided should not have a bond
            if p.type == self.part_type_to_be_glued:
                self.assertEqual(len(p.bonds), 0)

            for bond in p.bonds:
                # Bond type and partner type
                # part_type_after_glueing can have a bond to a vs or to a
                # non_virtual particle
                allowed_types = (self.part_type_after_glueing,
                                 self.part_type_to_attach_vs_to)
                self.assertIn(
                    p.type,
                    allowed_types,
                    msg=f"Particle {p.id} of type {p.type} should not have bonds, yet has {p.bonds}.")
                if p.type == self.part_type_after_glueing:
                    self.assertIn(bond[0], (self.bond_center, self.bond_vs))
                    # Bonds to virtual sites:
                    if bond[0] == self.bond_vs:
                        self.assertEqual(
                            system.part.by_id(bond[1]).type,
                            self.part_type_vs)
                    else:
                        self.assertEqual(
                            system.part.by_id(bond[1]).type,
                            self.part_type_to_attach_vs_to)
                elif p.type == self.part_type_to_attach_vs_to:
                    self.assertEqual(bond[0], self.bond_center)
                    self.assertEqual(
                        system.part.by_id(bond[1]).type,
                        self.part_type_after_glueing)

                # Collect bonds
                # Sort bond partners to make them unique independently of
                # which particle got the bond
                if bond[0] == self.bond_center:
                    bonds_centers.append(tuple(sorted([p.id, bond[1]])))
                else:
                    bonds_virtual.append(tuple(sorted([p.id, bond[1]])))

        # No duplicate bonds?
        self.assertEqual(len(bonds_centers), len(set(bonds_centers)))
        self.assertEqual(len(bonds_virtual), len(set(bonds_virtual)))

        # 1 bond between centers and one between vs and glued particle
        # per collision
        self.assertEqual(len(bonds_virtual), len(bonds_centers))

        # 1 virtual sites per bond?
        self.assertEqual(len(bonds_centers), len(virtual_sites))

        # no bonds on vs and vs particle type
        for p in virtual_sites:
            self.assertEqual(len(p.bonds), 0)
            self.assertEqual(p.type, self.part_type_vs)

        # Tidy
        system.non_bonded_inter[0, 0].lennard_jones.deactivate()

    def verify_triangle_binding(self, distance, first_bond, angle_res):
        system = self.system
        # Gather pairs
        n = len(system.part)
        angle_res = angle_res - 1

        expected_pairs = []
        for i in range(n):
            for j in range(i + 1, n, 1):
                if system.distance(system.part.by_id(i),
                                   system.part.by_id(j)) <= distance:
                    expected_pairs.append((i, j))

        # Find triangles
        # Each element is a particle id, a bond id and two bond partners in
        # ascending order
        expected_angle_bonds = []
        for i in range(n):
            for j in range(i + 1, n, 1):
                for k in range(j + 1, n, 1):
                    # Ref to particles
                    p_i, p_j, p_k = system.part.by_ids([i, j, k])

                    # Normalized distance vectors
                    d_ij = np.copy(p_j.pos - p_i.pos)
                    d_ik = np.copy(p_k.pos - p_i.pos)
                    d_jk = np.copy(p_k.pos - p_j.pos)
                    d_ij /= np.linalg.norm(d_ij)
                    d_ik /= np.linalg.norm(d_ik)
                    d_jk /= np.linalg.norm(d_jk)

                    if system.distance(p_i, p_j) <= distance and system.distance(
                            p_i, p_k) <= distance:
                        id_i = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(d_ij, d_ik)) * angle_res / np.pi))
                        expected_angle_bonds.append((i, id_i, j, k))

                    if system.distance(p_i, p_j) <= distance and system.distance(
                            p_j, p_k) <= distance:
                        id_j = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(-d_ij, d_jk)) * angle_res / np.pi))
                        expected_angle_bonds.append((j, id_j, i, k))
                    if system.distance(p_i, p_k) <= distance and system.distance(
                            p_j, p_k) <= distance:
                        id_k = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(-d_ik, -d_jk)) * angle_res / np.pi))
                        expected_angle_bonds.append((k, id_k, i, j))

        # Gather actual pairs and actual triangles
        found_pairs = []
        found_angle_bonds = []
        for i in range(n):
            for b in system.part.by_id(i).bonds:
                self.assertIn(
                    len(b), (2, 3), msg="There should only be 2- and 3-particle bonds")
                if len(b) == 2:
                    self.assertEqual(b[0]._bond_id, self.bond_center._bond_id)
                    found_pairs.append(tuple(sorted((i, b[1]))))
                elif len(b) == 3:
                    partners = sorted(b[1:])
                    found_angle_bonds.append(
                        (i, b[0]._bond_id, partners[0], partners[1]))

        # The order between expected and found bonds does not always match
        # because collisions occur in random order. Sort stuff
        found_pairs = sorted(found_pairs)
        found_angle_bonds = sorted(found_angle_bonds)
        expected_angle_bonds = sorted(expected_angle_bonds)
        self.assertEqual(found_pairs, expected_pairs)
        self.assertEqual(found_angle_bonds, expected_angle_bonds)


if __name__ == "__main__":
    ut.main()
