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
from espressomd.interactions import HarmonicBond, AngleHarmonic
import numpy as np
from random import shuffle


@utx.skipIfMissingFeatures("COLLISION_DETECTION")
class CollisionDetection(ut.TestCase):

    """Tests interface and functionality of the collision detection / dynamic binding"""

    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(seed=42)
    if espressomd.has_features("VIRTUAL_SITES_RELATIVE"):
        from espressomd.virtual_sites import VirtualSitesRelative
        s.virtual_sites = VirtualSitesRelative()

    H = HarmonicBond(k=5000, r_0=0.1)
    H2 = HarmonicBond(k=25000, r_0=0.02)
    s.bonded_inter.add(H)
    s.bonded_inter.add(H2)
    time_step = 0.001
    s.time_step = time_step
    s.cell_system.skin = 0.05
    s.min_global_cut = 0.112

    part_type_to_attach_vs_to = 0
    part_type_vs = 1
    part_type_to_be_glued = 2
    part_type_after_glueing = 3
    other_type = 5

    def get_state_set_state_consistency(self):
        state = self.s.collision_detection.get_params()
        self.s.collision_detection.set_params(**state)
        self.assertEqual(state, self.s.collision_detection.get_params())

    def test_00_interface_and_defaults(self):
        # Is it off by default
        self.assertEqual(self.s.collision_detection.mode, "off")
        # Make sure params cannot be set individually
        with self.assertRaises(Exception):
            self.s.collision_detection.mode = "bind_centers"

        # Verify exception throwing for unknown collision modes
        with self.assertRaises(Exception):
            self.s.collision_detection.set_params(mode=0)
            self.s.collision_detection.set_params(mode="blahblah")

        # That should work
        self.s.collision_detection.set_params(mode="off")
        self.assertEqual(self.s.collision_detection.mode, "off")

    def test_bind_centers(self):
        # Check that it leaves particles alone, when off
        self.s.collision_detection.set_params(mode="off")

        self.s.part.clear()
        self.s.part.add(pos=(0, 0, 0), id=0)
        self.s.part.add(pos=(0.1, 0, 0), id=1)
        self.s.part.add(pos=(0.1, 0.3, 0), id=2)
        self.s.integrator.run(1)
        self.assertEqual(self.s.part[0].bonds, ())
        self.assertEqual(self.s.part[1].bonds, ())
        self.assertEqual(self.s.part[2].bonds, ())

        # Check that it cannot be activated
        self.s.collision_detection.set_params(
            mode="bind_centers", distance=0.11, bond_centers=self.H)
        self.get_state_set_state_consistency()
        self.s.integrator.run(1, recalc_forces=True)
        bond0 = ((self.s.bonded_inter[0], 1),)
        bond1 = ((self.s.bonded_inter[0], 0),)
        self.assertTrue(
            self.s.part[0].bonds == bond0 or self.s.part[1].bonds == bond1)
        self.assertEqual(self.s.part[2].bonds, ())

        # Check that no additional bonds appear
        self.s.integrator.run(1)
        self.assertTrue(
            self.s.part[0].bonds == bond0 or self.s.part[1].bonds == bond1)
        self.assertEqual(self.s.part[2].bonds, ())

        # Check turning it off
        self.s.collision_detection.set_params(mode="off")
        self.get_state_set_state_consistency()
        self.assertEqual(self.s.collision_detection.mode, "off")

    def run_test_bind_at_point_of_collision_for_pos(self, *positions):
        positions = list(positions)
        shuffle(positions)
        self.s.part.clear()
        # Place particle which should not take part in collisions
        p = self.s.part.add(pos=(0.1, 0.3, 0))
        for pos in positions:
            p1 = self.s.part.add(pos=pos + (0, 0, 0))
            p2 = self.s.part.add(pos=pos + (0.1, 0, 0))
            if self.s.distance(p1, p) < 0.12 or self.s.distance(p2, p) < 0.12:
                raise Exception(
                    "Test particle too close to particle, which should not take part in collision")

        # 2 non-virtual + 2 virtual + one that doesn't take part
        expected_np = 4 * len(positions) + 1

        self.s.collision_detection.set_params(
            mode="bind_at_point_of_collision", distance=0.11, bond_centers=self.H, bond_vs=self.H2, part_type_vs=1, vs_placement=0.4)
        self.get_state_set_state_consistency()
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Integrate again and check that nothing has changed
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_state_after_bind_at_poc(expected_np)

        # Check that nothing explodes when the particles are moved.
        # In particular for parallel simulations
        self.s.thermostat.set_langevin(kT=0, gamma=0.01, seed=42)
        self.s.part[:].v = [0.05, 0.01, 0.15]
        self.s.integrator.run(3000)
        self.verify_state_after_bind_at_poc(expected_np)

    def verify_state_after_bind_at_poc(self, expected_np):
        self.assertEqual(len(self.s.part), expected_np)

        # At the end of test, this list should be empty
        parts_not_accounted_for = list(range(expected_np))

        # We traverse particles. We look for a vs with a bond to find the other vs.
        # From the two vs we find the two non-virtual particles
        for p in self.s.part:
            # Skip non-virtual
            if not p.virtual:
                continue
            # Skip vs that doesn't have a bond
            if p.bonds == ():
                continue
            # Parse the bond
            self.assertEqual(len(p.bonds), 1)
            # Bond type
            self.assertEqual(p.bonds[0][0], self.H2)
            # get partner
            p2 = self.s.part[p.bonds[0][1]]
            # Is that really a vs
            self.assertTrue(p2.virtual)
            # Get base particles
            base_p1 = self.s.part[p.vs_relative[0]]
            base_p2 = self.s.part[p2.vs_relative[0]]
            # Take note of accounted-for particles
            for _p in p, p2, base_p1, base_p2:
                parts_not_accounted_for.remove(_p.id)
            self.verify_bind_at_poc_pair(base_p1, base_p2, p, p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for), 1)
        p = self.s.part[parts_not_accounted_for[0]]
        self.assertFalse(p.virtual)
        self.assertEqual(p.bonds, ())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for, [])

    def verify_bind_at_poc_pair(self, p1, p2, vs1, vs2):
        bond_p1 = ((self.s.bonded_inter[0], p2.id),)
        bond_p2 = ((self.s.bonded_inter[0], p1.id),)
        self.assertTrue(p1.bonds == bond_p1 or p2.bonds == bond_p2)

        # Check for presence of vs
        # Check for bond between vs
        bond_vs1 = ((self.s.bonded_inter[1], vs2.id),)
        bond_vs2 = ((self.s.bonded_inter[1], vs1.id),)
        self.assertTrue(vs1.bonds == bond_vs1 or vs2.bonds == bond_vs2)

        # Vs properties
        self.assertTrue(vs1.virtual)
        self.assertTrue(vs2.virtual)

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
            expected_pos = self.s.part[rel_to].pos_folded + \
                self.s.collision_detection.vs_placement * dist_centers
            dist = expected_pos - p.pos_folded
            dist -= np.round(dist / self.s.box_l) * self.s.box_l
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

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "VIRTUAL_SITES_RELATIVE"])
    def test_bind_at_point_of_collision_random(self):
        """Integrate lj liquid and check that no double bonds are formed
           and the number of bonds fits the number of virtual sites

        """
        self.s.part.clear()

        # Add randomly placed particles
        self.s.part.add(pos=np.random.random((200, 3)))

        # Setup Lennard-Jones
        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=0.1, cutoff=2**(1. / 6) * 0.1, shift="auto")

        # Remove overlap between particles
        self.s.integrator.set_steepest_descent(
            f_max=0,
            gamma=1,
            max_displacement=0.001)
        while self.s.analysis.energy()["total"] > len(self.s.part):
            self.s.integrator.run(10)

        # Collision detection
        self.s.collision_detection.set_params(
            mode="bind_at_point_of_collision",
            distance=0.11,
            bond_centers=self.H,
            bond_vs=self.H2,
            part_type_vs=1,
            vs_placement=0.4)
        self.get_state_set_state_consistency()

        # Integrate lj liquid
        self.s.integrator.set_vv()
        self.s.integrator.run(5000)

        # Analysis
        virtual_sites = self.s.part.select(virtual=True)
        non_virtual = self.s.part.select(virtual=False)

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
                [self.s.part[vs_pair[0]].vs_relative[0],
                 self.s.part[vs_pair[1]].vs_relative[0]]))

            # Is there a corresponding bond?
            self.assertIn(base_particles, bonds)

        # Tidy
        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0)

    def run_test_glue_to_surface_for_pos(self, *positions):
        positions = list(positions)
        shuffle(positions)
        self.s.part.clear()
        # Place particle which should not take part in collisions
        # In this case, it is skipped, because it is of the wrong type,
        # even if it is within range for a collision
        self.s.part.add(pos=positions[0], type=self.other_type)
        for pos in positions:
            # Since this is non-symmetric, we randomize order
            if np.random.random() > .5:
                self.s.part.add(
                    pos=pos + (0, 0, 0), type=self.part_type_to_attach_vs_to)
                self.s.part.add(
                    pos=pos + (0.1, 0, 0), type=self.part_type_to_be_glued)
            else:
                self.s.part.add(
                    pos=pos + (0.1, 0, 0), type=self.part_type_to_be_glued)
                self.s.part.add(
                    pos=pos + (0, 0, 0), type=self.part_type_to_attach_vs_to)

        # 2 non-virtual + 1 virtual + one that doesn't take part
        expected_np = 3 * len(positions) + 1

        self.s.collision_detection.set_params(
            mode="glue_to_surface", distance=0.11, distance_glued_particle_to_vs=0.02, bond_centers=self.H, bond_vs=self.H2, part_type_vs=self.part_type_vs, part_type_to_attach_vs_to=self.part_type_to_attach_vs_to, part_type_to_be_glued=self.part_type_to_be_glued, part_type_after_glueing=self.part_type_after_glueing)
        self.get_state_set_state_consistency()
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)

        # Integrate again and check that nothing has changed
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_state_after_glue_to_surface(expected_np)

        # Check that nothing explodes, when the particles are moved.
        # In particular for parallel simulations
        self.s.thermostat.set_langevin(kT=0, gamma=0.01, seed=42)
        self.s.part[:].v = [0.05, 0.01, 0.15]
        self.s.integrator.run(3000)
        self.verify_state_after_glue_to_surface(expected_np)

    def verify_state_after_glue_to_surface(self, expected_np):
        self.assertEqual(len(self.s.part), expected_np)

        # At the end of test, this list should be empty
        parts_not_accounted_for = list(range(expected_np))

        # We traverse particles. We look for a vs, get base particle from there
        # and partner particle via bonds
        for p in self.s.part:
            # Skip non-virtual
            if not p.virtual:
                continue
            # The vs shouldn't have bonds
            self.assertEqual(p.bonds, ())

            # Get base particles
            base_p = self.s.part[p.vs_relative[0]]

            # Get bound particle
            # There is a bond between the base particle and the bound particle
            # but we have no guarantee, on where it is stored
            # 1. On the base particle of the vs
            p2 = None
            if len(base_p.bonds) == 1:
                self.assertEqual(base_p.bonds[0][0], self.H)
                p2 = self.s.part[base_p.bonds[0][1]]
            else:
                # We need to go through all particles to find it
                for candidate in self.s.part:
                    if candidate.id not in parts_not_accounted_for:
                        continue
                    if len(candidate.bonds) >= 1:
                        for b in candidate.bonds:
                            if b[0] == self.H and b[1] == base_p.id:
                                p2 = candidate
                if p2 is None:
                    raise Exception("Bound particle not found")
            # Take note of accounted-for particles
            parts_not_accounted_for.remove(base_p.id)
            parts_not_accounted_for.remove(p.id)
            parts_not_accounted_for.remove(p2.id)
            self.verify_glue_to_surface_pair(base_p, p, p2)
        # Check particle that did not take part in collision.
        self.assertEqual(len(parts_not_accounted_for), 1)
        p = self.s.part[parts_not_accounted_for[0]]
        self.assertFalse(p.virtual)
        self.assertEqual(p.type, self.other_type)
        self.assertEqual(p.bonds, ())
        parts_not_accounted_for.remove(p.id)
        self.assertEqual(parts_not_accounted_for, [])

    def verify_glue_to_surface_pair(self, base_p, vs, bound_p):
        # Check all types
        self.assertEqual(base_p.type, self.part_type_to_attach_vs_to)
        self.assertEqual(vs.type, self.part_type_vs)
        self.assertEqual(bound_p.type, self.part_type_after_glueing)

        # Bound particle should have a bond to vs. It can additionally have a bond
        # to the base particle
        bond_to_vs_found = 0
        for b in bound_p.bonds:
            if b[0] == self.H2:
                # bond to vs
                self.assertEqual(b, (self.H2, vs.id))
                bond_to_vs_found += 1
        self.assertEqual(bond_to_vs_found, 1)
        # Vs should not have a bond
        self.assertEqual(vs.bonds, ())

        # Vs properties
        self.assertTrue(vs.virtual)
        self.assertEqual(vs.vs_relative[0], base_p.id)

        # Distance vs,bound_p
        self.assertAlmostEqual(self.s.distance(vs, bound_p), 0.02, places=3)
        self.assertAlmostEqual(self.s.distance(base_p, bound_p), 0.1, places=3)
        self.assertAlmostEqual(self.s.distance(base_p, vs), 0.08, places=3)

        # base_p,vs,bound_p on a line
        self.assertGreater(np.dot(self.s.distance_vec(base_p, vs), self.s.distance_vec(base_p, bound_p))
                           / self.s.distance(base_p, vs) / self.s.distance(base_p, bound_p), 0.99)

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

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_glue_to_surface_random(self):
        """Integrate lj liquid and check that no double bonds are formed
           and the number of bonds fits the number of virtual sites

        """
        self.s.part.clear()

        # Add randomly placed particles
        self.s.part.add(pos=np.random.random((100, 3)),
                        type=100 * [self.part_type_to_attach_vs_to])
        self.s.part.add(pos=np.random.random(
            (100, 3)), type=100 * [self.part_type_to_be_glued])
        self.s.part.add(pos=np.random.random(
            (100, 3)), type=100 * [self.other_type])

        # Setup Lennard-Jones
        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=0.1, cutoff=2**(1. / 6) * 0.1, shift="auto")

        # Remove overlap between particles
        self.s.integrator.set_steepest_descent(
            f_max=0,
            gamma=1,
            max_displacement=0.001)
        while self.s.analysis.energy()["total"] > len(self.s.part):
            self.s.integrator.run(10)

        # Collision detection
        self.s.collision_detection.set_params(
            mode="glue_to_surface", distance=0.11, distance_glued_particle_to_vs=0.02, bond_centers=self.H, bond_vs=self.H2, part_type_vs=self.part_type_vs, part_type_to_attach_vs_to=self.part_type_to_attach_vs_to, part_type_to_be_glued=self.part_type_to_be_glued, part_type_after_glueing=self.part_type_after_glueing)
        self.get_state_set_state_consistency()

        # Integrate lj liquid
        self.s.integrator.set_vv()
        self.s.integrator.run(500)

        # Analysis
        virtual_sites = self.s.part.select(virtual=True)
        non_virtual = self.s.part.select(virtual=False)
        after_glueing = self.s.part.select(type=self.part_type_after_glueing)

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
                if p.type == self.part_type_after_glueing:
                    self.assertIn(bond[0], (self.H, self.H2))
                    # Bonds to virtual sites:
                    if bond[0] == self.H2:
                        self.assertEqual(
                            self.s.part[bond[1]].type,
                            self.part_type_vs)
                    else:
                        self.assertEqual(
                            self.s.part[bond[1]].type,
                            self.part_type_to_attach_vs_to)
                elif p.type == self.part_type_to_attach_vs_to:
                    self.assertEqual(bond[0], self.H)
                    self.assertEqual(
                        self.s.part[bond[1]].type,
                        self.part_type_after_glueing)
                else:
                    print(p.id, p.type, p.bonds)
                    raise Exception("Particle should not have bonds. ")

                # Collect bonds
                # Sort bond partners to make them unique independently of
                # which particle got the bond
                if bond[0] == self.H:
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
        self.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0, sigma=0, cutoff=0)

    def test_bind_three_particles(self):
        # Setup particles
        self.s.part.clear()
        dx = np.array((1, 0, 0))
        dy = np.array((0, 1, 0))
        a = np.array((0.499, 0.499, 0.499))
        b = a + 0.1 * dx
        c = a + 0.03 * dx + 0.03 * dy
        d = a + 0.03 * dx - 0.03 * dy
        e = a - 0.1 * dx

        self.s.part.add(id=0, pos=a)
        self.s.part.add(id=1, pos=b)
        self.s.part.add(id=2, pos=c)
        self.s.part.add(id=3, pos=d)
        self.s.part.add(id=4, pos=e)

        # Setup bonds
        res = 181
        for i in range(0, res, 1):
            self.s.bonded_inter[i + 2] = AngleHarmonic(
                bend=1, phi0=float(i) / (res - 1) * np.pi)
        cutoff = 0.11
        self.s.collision_detection.set_params(
            mode="bind_three_particles", bond_centers=self.H,
            bond_three_particles=2, three_particle_binding_angle_resolution=res, distance=cutoff)
        self.get_state_set_state_consistency()

        self.s.time_step = 1E-6
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_triangle_binding(cutoff, self.s.bonded_inter[2], res)
        # Make sure no extra bonds appear
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_triangle_binding(cutoff, self.s.bonded_inter[2], res)

        # Place the particles in two steps and make sure, the bonds are the
        # same
        self.s.part.clear()
        self.s.part.add(id=0, pos=a)
        self.s.part.add(id=2, pos=c)
        self.s.part.add(id=3, pos=d)
        self.s.integrator.run(1, recalc_forces=True)

        self.s.part.add(id=4, pos=e)
        self.s.part.add(id=1, pos=b)
        self.s.cell_system.set_domain_decomposition()
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_triangle_binding(cutoff, self.s.bonded_inter[2], res)
        self.s.cell_system.set_n_square()
        self.s.part[:].bonds = ()
        self.s.integrator.run(1, recalc_forces=True)
        self.verify_triangle_binding(cutoff, self.s.bonded_inter[2], res)
        self.s.time_step = self.time_step

    def verify_triangle_binding(self, distance, first_bond, angle_res):
        # Gather pairs
        n = len(self.s.part)
        angle_res = angle_res - 1

        expected_pairs = []
        for i in range(n):
            for j in range(i + 1, n, 1):
                if self.s.distance(self.s.part[i], self.s.part[j]) <= distance:
                    expected_pairs.append((i, j))

        # Find triangles
        # Each element is a particle id, a bond id and two bond partners in
        # ascending order
        expected_angle_bonds = []
        for i in range(n):
            for j in range(i + 1, n, 1):
                for k in range(j + 1, n, 1):
                    # Ref to particles
                    p_i = self.s.part[i]
                    p_j = self.s.part[j]
                    p_k = self.s.part[k]

                    # Normalized distance vectors
                    d_ij = np.copy(p_j.pos - p_i.pos)
                    d_ik = np.copy(p_k.pos - p_i.pos)
                    d_jk = np.copy(p_k.pos - p_j.pos)
                    d_ij /= np.linalg.norm(d_ij)
                    d_ik /= np.linalg.norm(d_ik)
                    d_jk /= np.linalg.norm(d_jk)

                    if self.s.distance(p_i, p_j) <= distance and self.s.distance(
                            p_i, p_k) <= distance:
                        id_i = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(d_ij, d_ik)) * angle_res / np.pi))
                        expected_angle_bonds.append((i, id_i, j, k))

                    if self.s.distance(p_i, p_j) <= distance and self.s.distance(
                            p_j, p_k) <= distance:
                        id_j = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(-d_ij, d_jk)) * angle_res / np.pi))
                        expected_angle_bonds.append((j, id_j, i, k))
                    if self.s.distance(p_i, p_k) <= distance and self.s.distance(
                            p_j, p_k) <= distance:
                        id_k = first_bond._bond_id + \
                            int(np.round(
                                np.arccos(np.dot(-d_ik, -d_jk)) * angle_res / np.pi))
                        expected_angle_bonds.append((k, id_k, i, j))

        # Gather actual pairs and actual triangles
        found_pairs = []
        found_angle_bonds = []
        for i in range(n):
            for b in self.s.part[i].bonds:
                if len(b) == 2:
                    self.assertEqual(b[0]._bond_id, self.H._bond_id)
                    found_pairs.append(tuple(sorted((i, b[1]))))
                elif len(b) == 3:
                    partners = sorted(b[1:])
                    found_angle_bonds.append(
                        (i, b[0]._bond_id, partners[0], partners[1]))
                else:
                    raise Exception(
                        "There should be only 2 and three particle bonds")

        # The order between expected and found bonds does not always match
        # because collisions occur in random order. Sort stuff
        found_pairs = sorted(found_pairs)
        found_angle_bonds = sorted(found_angle_bonds)
        expected_angle_bonds = sorted(expected_angle_bonds)
        self.assertEqual(expected_pairs, found_pairs)

        if not expected_angle_bonds == found_angle_bonds:
            # Verbose info
            print("expected:", expected_angle_bonds)
            missing = []
            for b in expected_angle_bonds:
                if b in found_angle_bonds:
                    found_angle_bonds.remove(b)
                else:
                    missing.append(b)
            print("missing", missing)
            print("extra:", found_angle_bonds)
            print()

        self.assertEqual(expected_angle_bonds, found_angle_bonds)

    def test_zz_serialization(self):
        self.s.collision_detection.set_params(
            mode="bind_centers", distance=0.11, bond_centers=self.H)
        reduce = self.s.collision_detection.__reduce__()
        res = reduce[0](reduce[1][0])
        self.assertEqual(res.__class__.__name__, "CollisionDetection")
        self.assertEqual(res.mode, "bind_centers")
        self.assertAlmostEqual(res.distance, 0.11, delta=1E-12)
        self.assertEqual(res.bond_centers, self.H)


if __name__ == "__main__":
    ut.main()
