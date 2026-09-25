#
# Copyright (C) 2013-2026 The ESPResSo project
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

import espressomd
import espressomd.reaction_methods

import numpy as np
import unittest as ut


class TestSingleReaction(ut.TestCase):

    def test_single_reaction(self):
        """
        Check a simple chemical reaction, the Monte Carlo acceptance rate
        and the configurational move probability for a given system state.
        """
        # create a reaction A -> 3 B + 4 C
        type_A = 0
        type_B = 1
        type_C = 2
        reaction = espressomd.reaction_methods.SingleReaction(
            gamma=2., reactant_types=[type_A], reactant_coefficients=[1],
            product_types=[type_B, type_C], product_coefficients=[3, 4])

        self.assertEqual(reaction.nu_bar, 6)
        self.assertEqual(reaction.trial_moves, 0)
        self.assertEqual(reaction.accepted_moves, 0)

        # check acceptance rate
        for trial_moves in range(1, 5):
            for accepted_moves in range(0, 5):
                reaction.trial_moves = trial_moves
                reaction.accepted_moves = accepted_moves
                ref_rate = float(accepted_moves) / float(trial_moves)
                self.assertAlmostEqual(
                    reaction.get_acceptance_rate(), ref_rate, delta=1e-12)

        # check factorial expression
        Method = espressomd.reaction_methods.ReactionAlgorithm
        ln_fact = Method._ln_factorial_Ni0_div_factorial_Ni0_nu_i
        calc_factorial_expression = espressomd.reaction_methods.ReactionEnsemble.calculate_factorial_expression
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    # system contains i x A, j x B, and k x C
                    p_numbers = {type_A: i, type_B: j, type_C: k}
                    val = calc_factorial_expression(reaction, p_numbers)
                    ref = ln_fact(i, -1) + ln_fact(j, 3) + ln_fact(k, 4)
                    self.assertAlmostEqual(val, ref, delta=1e-12)


class ReactionMethods(ut.TestCase):

    """Test the reaction methods interface."""

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    node_grid = np.copy(system.cell_system.node_grid)

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.cell_system.node_grid = self.node_grid

    def tearDown(self):
        self.system.part.clear()

    def check_interface(self, method, kT, exclusion_range,
                        gamma, exclusion_radius_per_type, search_algorithm):
        def check_reaction_parameters(reactions, parameters):
            for reaction, params in zip(reactions, parameters):
                for key in reaction.required_keys():
                    if isinstance(params[key], float):
                        self.assertAlmostEqual(
                            getattr(reaction, key), params[key], delta=1e-10)
                    else:
                        self.assertEqual(getattr(reaction, key), params[key])

        def count_by_type(types):
            return [self.system.number_of_particles(type=x) for x in types]

        # check acceptance rate
        for tried_moves in range(1, 5):
            for accepted_moves in range(0, 5):
                method.m_tried_configurational_MC_moves = tried_moves
                method.m_accepted_configurational_MC_moves = accepted_moves
                ref_rate = float(accepted_moves) / float(tried_moves)
                self.assertAlmostEqual(
                    method.get_acceptance_rate_configurational_moves(),
                    ref_rate, delta=1e-10
                )
        method.m_tried_configurational_MC_moves = 0
        method.m_accepted_configurational_MC_moves = 0

        # add reactions
        reaction_forward = {
            "gamma": gamma,
            "reactant_types": [5],
            "reactant_coefficients": [1],
            "product_types": [2, 3],
            "product_coefficients": [1, 1],
            "default_charges": {5: 0, 2: 0, 3: 0},
        }
        reaction_backward = {
            "gamma": 1. / gamma,
            "reactant_types": reaction_forward["product_types"],
            "reactant_coefficients": reaction_forward["product_coefficients"],
            "product_types": reaction_forward["reactant_types"],
            "product_coefficients": reaction_forward["reactant_coefficients"],
            "default_charges": reaction_forward["default_charges"],
        }

        if isinstance(method, espressomd.reaction_methods.ConstantpHEnsemble):
            method.add_reaction(gamma=reaction_forward["gamma"],
                                reactant_types=reaction_forward["reactant_types"],
                                product_types=reaction_forward["product_types"],
                                default_charges=reaction_forward["default_charges"])
        else:
            method.add_reaction(**reaction_forward)
        reaction_parameters = (reaction_forward, reaction_backward)

        # check getters and setters
        self.assertAlmostEqual(method.kT, kT, delta=1e-10)
        self.assertAlmostEqual(
            method.exclusion_range,
            exclusion_range,
            delta=1e-10)
        self.assertEqual(method.search_algorithm, search_algorithm)
        if not isinstance(method, espressomd.reaction_methods.WidomInsertion):
            self.assertEqual(
                list(method.exclusion_radius_per_type.keys()), [1])
            self.assertAlmostEqual(
                method.exclusion_radius_per_type[1],
                exclusion_radius_per_type[1],
                delta=1e-10)
            method.exclusion_radius_per_type = {2: 0.2}
            self.assertEqual(
                list(method.exclusion_radius_per_type.keys()), [2])
            self.assertAlmostEqual(
                method.exclusion_radius_per_type[2], 0.2, delta=1e-10)
            exclusion_radius_per_type = {2: 0.2}
            method.exclusion_range = 0.8
            self.assertAlmostEqual(method.exclusion_range, 0.8, delta=1e-10)
            method.exclusion_range = 0.8
        self.assertAlmostEqual(
            method.get_volume(), self.system.volume(), delta=1e-10)
        self.assertEqual(method.get_non_interacting_type(), 100)
        method.set_non_interacting_type(type=9)
        self.assertEqual(method.get_non_interacting_type(), 9)
        if isinstance(method, espressomd.reaction_methods.ConstantpHEnsemble):
            self.assertAlmostEqual(method.constant_pH, 10., delta=1e-10)
            method.constant_pH = 8.
            self.assertAlmostEqual(method.constant_pH, 8., delta=1e-10)
        self.assertFalse(method.displacement_mc_move_for_particles_of_type(
            type_mc=0, particle_number_to_be_changed=0))
        self.assertFalse(method.displacement_mc_move_for_particles_of_type(
            type_mc=0, particle_number_to_be_changed=100000))
        method.particle_inside_exclusion_range_touched = True
        self.assertTrue(method.particle_inside_exclusion_range_touched)
        method.particle_inside_exclusion_range_touched = False
        self.assertFalse(method.particle_inside_exclusion_range_touched)

        # check constraints
        method.set_wall_constraints_in_z_direction(
            slab_start_z=0.1, slab_end_z=0.9)
        offsets = method.get_wall_constraints_in_z_direction()
        self.assertAlmostEqual(offsets[0], 0.1, delta=1e-10)
        self.assertAlmostEqual(offsets[1], 0.9, delta=1e-10)
        method.remove_constraint()

        # check status
        status = method.get_status()
        self.assertEqual(status["kT"], kT)
        self.assertEqual(status["exclusion_range"], exclusion_range)
        self.assertEqual(
            status["exclusion_radius_per_type"],
            exclusion_radius_per_type)
        self.assertEqual(len(status["reactions"]), 2)
        for reaction_flat, params in zip(
                status["reactions"], reaction_parameters):
            for key in reaction_flat:
                if isinstance(params[key], float):
                    self.assertAlmostEqual(
                        reaction_flat[key], params[key], delta=1e-10)
                else:
                    self.assertEqual(reaction_flat[key], params[key])

        # check reactions
        reactions = method.reactions
        self.assertEqual(len(reactions), 2)
        check_reaction_parameters(method.reactions, reaction_parameters)

        # check reactions after unsuccessful parameter change
        with self.assertRaises(ValueError):
            method.change_reaction_constant(reaction_id=0, gamma=0.)
        check_reaction_parameters(method.reactions, reaction_parameters)

        # check reactions after successful parameter change
        new_gamma = 634.
        reaction_forward["gamma"] = new_gamma
        reaction_backward["gamma"] = 1. / new_gamma
        method.change_reaction_constant(reaction_id=0, gamma=new_gamma)
        check_reaction_parameters(method.reactions, reaction_parameters)
        status = method.get_status()
        self.assertAlmostEqual(
            status["reactions"][0]["gamma"],
            reaction_forward["gamma"],
            delta=1e-10)
        self.assertAlmostEqual(
            status["reactions"][1]["gamma"],
            reaction_backward["gamma"],
            delta=1e-10)

        # check particle deletion on a worker node
        self.system.part.clear()
        p1, _, p3 = self.system.part.add(
            pos=3 * [(-1., -1., -1.)], type=[5, 2, 3])
        if isinstance(method, espressomd.reaction_methods.WidomInsertion):
            potential_energy = method.calculate_particle_insertion_potential_energy(
                reaction_id=0)
            self.assertEqual(potential_energy, 0.)
        if self.system.cell_system.get_state()["n_nodes"] > 1:
            assert set(self.system.part.all().node) != {0}
        self.assertEqual(count_by_type([5, 2, 3, 0]), [1, 1, 1, 0])
        method.delete_particle(p_id=p3.id)
        self.assertEqual(count_by_type([5, 2, 3, 0]), [1, 1, 0, 0])
        self.assertEqual(len(self.system.part), 2)
        p1.remove()
        self.assertEqual(count_by_type([5, 2, 3, 0]), [0, 1, 0, 0])
        self.assertEqual(len(self.system.part), 1)
        self.system.part.clear()
        self.assertEqual(count_by_type([5, 2, 3, 0]), [0, 0, 0, 0])

        # check particle hiding functions
        self.system.part.clear()
        p1, p2, p3 = self.system.part.add(pos=3 * [(0., 0., 0.)], type=3 * [5])
        method._hide_particle(p1.id)
        self.assertEqual(p1.type, method.non_interacting_type)
        self.assertEqual(p2.type, 5)
        self.assertEqual(p3.type, 5)
        method._hide_particles(pids=[p2.id, p3.id], ptype=5)
        self.assertEqual(p1.type, method.non_interacting_type)
        self.assertEqual(p2.type, method.non_interacting_type)
        self.assertEqual(p3.type, method.non_interacting_type)

    def test_reaction_interface(self):
        params = {"exclusion_range": 0.8,
                  "exclusion_radius_per_type": {1: 0.1}}

        with self.subTest(msg="reaction ensemble"):
            method = espressomd.reaction_methods.ReactionEnsemble(
                kT=1.4, seed=12, search_algorithm="order_n", system=self.system, **params)
            self.check_interface(method, kT=1.4, gamma=1.2,
                                 search_algorithm="order_n", **params)

        if espressomd.has_features("ELECTROSTATICS"):
            with self.subTest(msg="constant pH ensemble"):
                method = espressomd.reaction_methods.ConstantpHEnsemble(
                    kT=1.5, seed=14, search_algorithm="parallel", constant_pH=10., system=self.system,
                    **params)
                self.check_interface(method, kT=1.5, gamma=1.2,
                                     search_algorithm="parallel", **params)

        with self.subTest(msg="Widom insertion"):
            method = espressomd.reaction_methods.WidomInsertion(
                kT=1.6, seed=16, system=self.system)
            self.check_interface(method, kT=1.6, gamma=1., exclusion_range=0.,
                                 exclusion_radius_per_type={},
                                 search_algorithm=None)

    def test_displacement_mc_move(self):
        ref_pos = [(0.1, 0.2, 0.3), (0.4, 0.5, 0.6)]
        ref_vel = [(10.0, 20.0, 30.0), (-10.0, -20.0, -30.0)]
        self.system.box_l = [1., 1., 1.]
        self.system.part.add(id=[0, 1], pos=ref_pos, v=ref_vel)

        r_algo = espressomd.reaction_methods.ReactionEnsemble(
            seed=42, kT=1., exclusion_range=0., system=self.system)
        r_algo.exclusion.exclusion_range = 1.
        self.assertFalse(r_algo.exclusion_range_touched)
        r_algo._displacement_mc_move(0, 2)
        self.assertTrue(r_algo.exclusion_range_touched)

        self.assertEqual(len(r_algo._particle_changes["created"]), 0)
        self.assertEqual(len(r_algo._particle_changes["hidden"]), 0)
        for change in r_algo._particle_changes["changed"]:
            pid = change["pid"]
            self.assertIn(pid, (0, 1))
            ref_old_pos = ref_pos[pid]
            ref_old_vel = ref_vel[pid]
            p = self.system.part.by_id(pid)
            new_pos = np.copy(p.pos)
            new_vel = np.copy(p.v)
            np.testing.assert_allclose(np.copy(change["pos"]), ref_old_pos)
            np.testing.assert_allclose(np.copy(change["v"]), ref_old_vel)
            self.assertTrue(np.all(new_pos >= 0.))
            self.assertTrue(np.all(new_pos <= 1.))
            self.assertGreaterEqual(np.linalg.norm(new_pos - ref_old_pos), 0.1)
            self.assertGreaterEqual(np.linalg.norm(new_vel - ref_old_vel), 10.)

    def test_displacement_mc_move_for_particles_of_type(self):
        ref_pos = [(0.1, 0.2, 0.3), (0.6, 0.7, 0.8)]
        self.system.box_l = [1., 1., 1.]
        self.system.part.add(id=[0, 1], pos=ref_pos)

        type_A = 0
        type_B = 1
        type_C = 2
        r_algo = espressomd.reaction_methods.ReactionEnsemble(
            seed=42, kT=1., exclusion_range=0., system=self.system)

        def displacement_move(t, n):
            return r_algo.displacement_mc_move_for_particles_of_type(t, n)

        # check impossible MC moves
        self.assertFalse(displacement_move(type_C, 1))
        self.assertFalse(displacement_move(type_B, 2))
        self.assertFalse(displacement_move(type_A, 0))
        with self.assertRaises(ValueError):
            displacement_move(type_A, -2)
        with self.assertRaises(ValueError):
            displacement_move(-2, 1)

        # force all MC moves to be rejected by picking particles inside their
        # exclusion radius
        r_algo.exclusion.exclusion_range = 1.
        self.assertFalse(displacement_move(type_A, 2))
        # check none of the particles moved
        for pid in (0, 1):
            ref_old_pos = ref_pos[pid]
            p = self.system.part.by_id(pid)
            np.testing.assert_allclose(np.copy(p.pos), ref_old_pos)
        # check that bookkeeping was reset
        self.assertFalse(r_algo._particle_changes["created"])
        self.assertFalse(r_algo._particle_changes["changed"])
        self.assertFalse(r_algo._particle_changes["hidden"])

        # force a MC move to be accepted by using a constant Hamiltonian
        r_algo.exclusion.exclusion_range = 0.
        self.assertTrue(displacement_move(type_A, 1))
        distances = [0.0, 0.0]
        # check that only one particle moved
        for pid in (0, 1):
            p = self.system.part.by_id(pid)
            distances[pid] = np.linalg.norm(ref_pos[pid] - p.pos)
        self.assertLessEqual(min(distances[0], distances[1]), 1e-10)
        self.assertGreaterEqual(max(distances[0], distances[1]), 0.1)
        # check that bookkeeping was reset
        self.assertFalse(r_algo._particle_changes["created"])
        self.assertFalse(r_algo._particle_changes["changed"])
        self.assertFalse(r_algo._particle_changes["hidden"])

    def test_constraints(self):
        box_l = np.array([0.5, 0.4, 0.7])
        origin = np.zeros(3)
        self.system.box_l = box_l
        r_algo = espressomd.reaction_methods.ReactionEnsemble(
            seed=40, kT=1., exclusion_range=0., system=self.system)

        # cubic case
        positions = r_algo.get_random_positions_in_box(100)
        for pos in positions:
            self.assertTrue(np.all(pos <= box_l))
            self.assertTrue(np.all(pos >= origin))

        # slab case
        start_z, end_z = (0.2, 0.6)
        slab_lower = np.array([0.0, 0.0, start_z])
        slab_upper = np.array([box_l[0], box_l[1], end_z])
        r_algo.set_wall_constraints_in_z_direction(start_z, end_z)
        slab_params = r_algo.get_wall_constraints_in_z_direction()
        self.assertAlmostEqual(slab_params[0], start_z, delta=1e-10)
        self.assertAlmostEqual(slab_params[1], end_z, delta=1e-10)
        positions = r_algo.get_random_positions_in_box(100)
        for pos in positions:
            self.assertTrue(np.all(pos <= slab_upper))
            self.assertTrue(np.all(pos >= slab_lower))

        # cylindrical case
        cyl_x, cyl_y, radius = (0.2, 0.1, 0.2)
        r_algo.set_cylindrical_constraint_in_z_direction(cyl_x, cyl_y, radius)
        positions = r_algo.get_random_positions_in_box(100)
        for pos in positions:
            z = pos[2]
            r = np.linalg.norm([pos[0] - cyl_x, pos[1] - cyl_y])
            self.assertLessEqual(r, radius)
            self.assertLessEqual(z, box_l[2])
            self.assertGreaterEqual(z, 0.0)

    def test_exclusion_agorithm(self):
        type_A = 0
        type_B = 1
        self.system.box_l = [1., 1., 1.]
        _, p1 = self.system.part.add(pos=[(0.5, 0.5, 0.5), (0.7, 0.7, 0.7)],
                                     type=[type_A, type_B])
        r_algo = espressomd.reaction_methods.ReactionEnsemble(
            seed=40, kT=1., exclusion_range=0., system=self.system)

        # new positions will always be in the excluded range if the sum of the
        # radii of both particle types is larger than box length (radii take
        # precedence over the default exclusion range)
        r_algo.exclusion.exclusion_range = 0.
        r_algo.exclusion.exclusion_radius_per_type = {type_A: 0.1, type_B: 2.}
        r_algo.exclusion_range_touched = False
        r_algo._displacement_mc_move(type_B, 1)
        self.assertTrue(r_algo.exclusion_range_touched)
        # also check private implementations
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=p1.id, ptype=type_B)
        self.assertTrue(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=1)
        self.assertTrue(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range_any(pids=[1], ptype=type_B)
        self.assertTrue(r_algo.exclusion_range_touched)

        # new positions will never be in the excluded range if the exclusion
        # radius of either particle is 0
        r_algo.exclusion.exclusion_range = 0.
        r_algo.exclusion.exclusion_radius_per_type = {type_A: 0., type_B: 2.}
        r_algo.exclusion_range_touched = False
        r_algo._displacement_mc_move(type_B, 1)
        self.assertFalse(r_algo.exclusion_range_touched)
        # also check private implementations
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=p1.id, ptype=type_B)
        self.assertFalse(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=1)
        self.assertFalse(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range_any(pids=[1], ptype=type_B)
        self.assertFalse(r_algo.exclusion_range_touched)

        # new positions will never be accepted if the exclusion range is larger
        # than box length and particles don't define radii to override it
        r_algo.exclusion.exclusion_range = 2.
        r_algo.exclusion.exclusion_radius_per_type = {type_A: 0.}
        r_algo.exclusion_range_touched = False
        r_algo._displacement_mc_move(type_B, 1)
        self.assertTrue(r_algo.exclusion_range_touched)
        # also check private implementations
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=p1.id, ptype=type_B)
        self.assertTrue(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range(pid=1)
        self.assertTrue(r_algo.exclusion_range_touched)
        r_algo.exclusion_range_touched = False
        r_algo._check_exclusion_range_any(pids=[1], ptype=type_B)
        self.assertTrue(r_algo.exclusion_range_touched)

    def test_exceptions(self):
        self.system.part.add(pos=3 * [(0., 0., 0.)], id=[0, 2, 4])
        single_reaction_params = {
            "gamma": 1.,
            "reactant_types": [4],
            "reactant_coefficients": [1],
            "product_types": [2, 3],
            "product_coefficients": [1, 4],
        }
        reaction_params = {
            "default_charges": {2: 0, 3: 0, 4: 0},
            **single_reaction_params
        }
        widom = espressomd.reaction_methods.WidomInsertion(
            kT=1., seed=12, system=self.system)
        method = espressomd.reaction_methods.ReactionEnsemble(
            kT=1.5, exclusion_range=0.1, seed=12, system=self.system,
            exclusion_radius_per_type={1: 0.1})
        method.add_reaction(**reaction_params)
        widom.add_reaction(**reaction_params)

        # check invalid reactions
        err_msg = "number of types and coefficients have to match"
        with self.assertRaisesRegex(ValueError, f"reactants: {err_msg}"):
            method.add_reaction(**{**reaction_params, "reactant_types": []})
        with self.assertRaisesRegex(ValueError, f"products: {err_msg}"):
            method.add_reaction(**{**reaction_params, "product_types": []})
        with self.assertRaisesRegex(ValueError, "gamma"):
            method.add_reaction(**{**reaction_params, "gamma": 0.})
        with self.assertRaisesRegex(ValueError, "gamma"):
            method.add_reaction(**{**reaction_params, "gamma": -2.})

        # check charge conservation
        err_msg = "Reaction system is not charge neutral"
        with self.assertRaisesRegex(ValueError, err_msg):
            method.add_reaction(default_charges={2: 8, 3: 0, 4: -50},
                                **single_reaction_params)
        with self.assertRaisesRegex(ValueError, err_msg):
            method.add_reaction(default_charges={2: 1, 3: 0, 4: 1 + 1e-10},
                                **single_reaction_params)
        with self.assertRaisesRegex(TypeError, "needs to be a dict"):
            method.add_reaction(default_charges=(1, 2),
                                **single_reaction_params)

        # check rollback
        self.assertEqual(len(method.reactions), 2)

        # check invalid reaction id exceptions
        # (note: reaction index = 2 * reaction id)
        for i in [-2, -1, 1, 2, 3]:
            with self.assertRaisesRegex(IndexError, f"No reaction with id {i}"):
                method.delete_reaction(reaction_id=i)
            with self.assertRaisesRegex(IndexError, f"No reaction with id {2 * i}"):
                method.get_acceptance_rate_reaction(reaction_id=2 * i)
        with self.assertRaisesRegex(IndexError, "No reaction with id 1"):
            method.change_reaction_constant(reaction_id=1, gamma=1.)

        # check constraint exceptions
        set_cyl_constraint = method.set_cylindrical_constraint_in_z_direction
        set_slab_constraint = method.set_wall_constraints_in_z_direction
        get_slab_constraint = method.get_wall_constraints_in_z_direction
        err_msg = "no slab constraint is currently active"
        with self.assertRaisesRegex(RuntimeError, err_msg):
            get_slab_constraint()
        set_slab_constraint(slab_start_z=0.1, slab_end_z=0.9)
        method.remove_constraint()
        with self.assertRaisesRegex(RuntimeError, err_msg):
            get_slab_constraint()

        # check invalid constraints
        with self.assertRaisesRegex(ValueError, "center_x is outside the box"):
            set_cyl_constraint(center_x=100., center_y=1., radius=1.)
        with self.assertRaisesRegex(ValueError, "center_x is outside the box"):
            set_cyl_constraint(center_x=-10., center_y=1., radius=1.)
        with self.assertRaisesRegex(ValueError, "center_y is outside the box"):
            set_cyl_constraint(center_y=100., center_x=1., radius=1.)
        with self.assertRaisesRegex(ValueError, "center_y is outside the box"):
            set_cyl_constraint(center_y=-10., center_x=1., radius=1.)
        with self.assertRaisesRegex(ValueError, "radius is invalid"):
            set_cyl_constraint(center_x=1., center_y=1., radius=-1.)
        with self.assertRaisesRegex(ValueError, "slab_start_z is outside the box"):
            set_slab_constraint(slab_start_z=100., slab_end_z=1.)
        with self.assertRaisesRegex(ValueError, "slab_start_z is outside the box"):
            set_slab_constraint(slab_start_z=-10., slab_end_z=1.)
        with self.assertRaisesRegex(ValueError, "slab_end_z is outside the box"):
            set_slab_constraint(slab_end_z=100., slab_start_z=1.)
        with self.assertRaisesRegex(ValueError, "slab_end_z is outside the box"):
            set_slab_constraint(slab_end_z=-10., slab_start_z=1.)
        with self.assertRaisesRegex(ValueError, "slab_end_z must be >= slab_start_z"):
            set_slab_constraint(slab_start_z=10., slab_end_z=1.)

        # check exceptions for missing particles
        with self.assertRaisesRegex(RuntimeError, "Particle id is greater than the max seen particle id"):
            method.delete_particle(p_id=6)
        with self.assertRaisesRegex(RuntimeError, "Particle with id 3 not found"):
            method.delete_particle(p_id=3)
        with self.assertRaisesRegex(ValueError, "Invalid particle id: -2"):
            method.delete_particle(p_id=-2)
        with self.assertRaisesRegex(RuntimeError, "Trying to remove some non-existing particles from the system via the inverse Widom scheme"):
            widom.calculate_particle_insertion_potential_energy(reaction_id=0)
        with self.assertRaisesRegex(RuntimeError, "No search algorithm for WidomInsertion"):
            widom.search_algorithm = "order_n"
        with self.assertRaisesRegex(RuntimeError, "No search algorithm for WidomInsertion"):
            widom.exclusion_range = 1.
        with self.assertRaisesRegex(RuntimeError, "No search algorithm for WidomInsertion"):
            widom.exclusion_radius_per_type = {1: 2.}

        # check other exceptions
        err_msg = r"Only the following keys can be given as keyword arguments: \[.+\], got \[.+\] \(unknown \['x'\]\)"
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_methods.SingleReaction(
                x=1, **single_reaction_params)
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_methods.ReactionEnsemble(
                kT=1., exclusion_range=1., seed=12, x=1, system=self.system)
        with self.assertRaisesRegex(ValueError, r"\(missing ..reactant_coefficients..\)"):
            espressomd.reaction_methods.SingleReaction(
                **{k: v for k, v in single_reaction_params.items() if k != "reactant_coefficients"})
        with self.assertRaisesRegex(ValueError, r"\(unknown ..unknown..\)"):
            espressomd.reaction_methods.SingleReaction(
                **single_reaction_params, unknown=5)
        if espressomd.has_features("ELECTROSTATICS"):
            cph_exception_context = self.assertRaisesRegex(ValueError, err_msg)
        else:
            cph_exception_context = self.assertRaisesRegex(
                espressomd.code_features.FeaturesError,
                "Missing features ELECTROSTATICS")
        with cph_exception_context:
            espressomd.reaction_methods.ConstantpHEnsemble(
                kT=1., exclusion_range=1., seed=12, x=1, constant_pH=2, system=self.system)
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_methods.WidomInsertion(
                kT=1., seed=12, x=1, system=self.system)
        with self.assertRaisesRegex(ValueError, "Invalid value for 'kT'"):
            espressomd.reaction_methods.ReactionEnsemble(
                kT=-1., exclusion_range=1., seed=12, system=self.system)
        if not espressomd.has_features("ELECTROSTATICS"):
            with self.assertRaisesRegex(RuntimeError, "Charge assigned to type 4 is non-zero, but ELECTROSTATICS is not compiled in"):
                tmp_method = espressomd.reaction_methods.ReactionEnsemble(
                    kT=1., exclusion_range=1., seed=12, system=self.system)
                tmp_method.add_reaction(
                    default_charges={2: -1, 3: 0, 4: -1},
                    **single_reaction_params)
            self.assertEqual(len(tmp_method.reactions), 0)  # check rollback
        with self.assertRaisesRegex(ValueError, "Parameter 'particle_number_to_be_changed' must be >= 0"):
            method.displacement_mc_move_for_particles_of_type(
                type_mc=0, particle_number_to_be_changed=-1)
        with self.assertRaisesRegex(ValueError, "Parameter 'type_mc' must be >= 0"):
            method.displacement_mc_move_for_particles_of_type(
                type_mc=-1, particle_number_to_be_changed=1)
        with self.assertRaisesRegex(RuntimeError, "cannot be instantiated"):
            espressomd.reaction_methods.ReactionAlgorithm()
        with self.assertRaisesRegex(ValueError, "Invalid type: -1"):
            method.set_non_interacting_type(type=-1)
        with self.assertRaisesRegex(NotImplementedError, "Derived classes must implement this method"):
            super(type(method), method).calculate_log_acceptance_probability(
                method.reactions[0], 0., {})
        with self.assertRaisesRegex(RuntimeError, "Types may not be negative"):
            method.non_interacting_type = -2
            method.count_number_of_particles_per_type()
        method.non_interacting_type = 100
        with self.assertRaisesRegex(RuntimeError, "unknown method 'unknown'"):
            method._helper.call_method("unknown")
        with self.assertRaisesRegex(RuntimeError, "Reaction methods do not support checkpointing"):
            method._helper._serialize()
        with self.assertRaisesRegex(NotImplementedError, "was removed in release 5.1"):
            self.system.setup_type_map()

        # check invalid exclusion ranges and radii
        with self.assertRaisesRegex(ValueError, "Invalid value for exclusion range"):
            espressomd.reaction_methods.ReactionEnsemble(
                kT=1., seed=12, exclusion_range=-1., system=self.system)
        with self.assertRaisesRegex(ValueError, r"Invalid exclusion radius for type 1: radius -0\.1[0]*"):
            espressomd.reaction_methods.ReactionEnsemble(
                kT=1., seed=12, exclusion_range=1., exclusion_radius_per_type={1: -0.1}, system=self.system)
        with self.assertRaisesRegex(ValueError, "Unknown search algorithm 'unknown'"):
            espressomd.reaction_methods.ReactionEnsemble(
                kT=1., seed=12, exclusion_range=1., search_algorithm="unknown", system=self.system)
        method = espressomd.reaction_methods.ReactionEnsemble(
            kT=1., exclusion_range=1., seed=12, exclusion_radius_per_type={1: 0.1}, system=self.system)
        with self.assertRaisesRegex(ValueError, r"Invalid exclusion radius for type 2: radius -0\.1[0]*"):
            method.exclusion_radius_per_type = {2: -0.1}
        self.assertEqual(list(method.exclusion_radius_per_type.keys()), [1])

    def test_change_constant_of_second_reaction(self):
        """
        Check reaction indexing. Some methods use the reaction index, while
        others use reaction id (internally, reaction index = 2 * reaction id).
        """
        method = espressomd.reaction_methods.ReactionEnsemble(
            kT=1., exclusion_range=1., seed=12, system=self.system)

        gamma0 = 2.   # initial forward gamma of the 1st reaction
        gamma1 = 5.   # initial forward gamma of the 2nd reaction

        # 1st logical reaction -> m_reactions indices 0 (forward), 1 (backward)
        method.add_reaction(
            gamma=gamma0,
            reactant_types=[0], reactant_coefficients=[1],
            product_types=[1, 2], product_coefficients=[1, 1],
            default_charges={0: 0, 1: 0, 2: 0})

        # 2nd logical reaction -> m_reactions indices 2 (forward), 3 (backward)
        method.add_reaction(
            gamma=gamma1,
            reactant_types=[3], reactant_coefficients=[1],
            product_types=[4, 5], product_coefficients=[1, 1],
            default_charges={3: 0, 4: 0, 5: 0})

        reactions = method.get_status()["reactions"]
        self.assertEqual(len(reactions), 4)
        # sanity: the flat container is [F0, B0, F1, B1] as expected
        self.assertAlmostEqual(reactions[0]["gamma"], gamma0, delta=1e-10)
        self.assertAlmostEqual(reactions[2]["gamma"], gamma1, delta=1e-10)

        # Per the documented convention, the SECOND added reaction is addressed
        # with reaction_id=1.  This must succeed and update indices 2 and 3.
        new_gamma = 17.
        method.change_reaction_constant(reaction_id=1, gamma=new_gamma)

        reactions = method.get_status()["reactions"]
        # second reaction forward/backward updated ...
        self.assertAlmostEqual(reactions[2]["gamma"], new_gamma, delta=1e-10)
        self.assertAlmostEqual(
            reactions[3]["gamma"], 1. / new_gamma, delta=1e-10)
        # ... and the first reaction left untouched
        self.assertAlmostEqual(reactions[0]["gamma"], gamma0, delta=1e-10)
        self.assertAlmostEqual(reactions[1]["gamma"], 1. / gamma0, delta=1e-10)

    @ut.skipIf(np.prod(node_grid) != 2, "needs 2 MPI ranks")
    def test_search_algorithm_order_n(self):
        """
        Test exclusion range larger than ghost layer. The order N range check
        must detect a candidate particle that is owned by a different MPI rank
        and is not a ghost on the inserted particle's rank. On the buggy code
        the distance loop runs only on the inserted particle's owning rank and
        resolves the candidate via get_local_particle (local + ghost only),
        so a remote, non-ghost particle within exclusion_range is silently
        missed and the flag stays false. The correct behavior is that the
        overlap is detected on any number of ranks.
        """
        # box 1x1x1; with 2 ranks the node grid is {2, 1, 1}, splitting the box
        # along x. The ghost layer (max_range) is then about one cell wide in x,
        # so a particle sitting 0.05 from an x-boundary is real on exactly one
        # rank and NOT a ghost on the other rank.
        self.system.box_l = [1., 1., 1.]

        # particle 0 -> owned by the rank holding small x, particle 1 -> owned
        # by the rank holding large x; neither is a ghost on the other rank.
        self.system.part.add(id=0, pos=[0.05, 0.5, 0.5])
        self.system.part.add(id=1, pos=[0.95, 0.5, 0.5])

        # precondition: the ghost layer is much narrower than the distance of
        # either particle to the far domain, so the particles are genuinely not
        # mutual ghosts (the condition that makes the bug observable)
        assert self.system.cell_system.call_method("get_max_range")[0] < 0.45

        # exclusion_range (0.5) larger than the minimum-image distance between
        # the two particles (mi-distance is 0.1) so an exclusion-range violation
        # exists, but the candidate is on a remote rank beyond the ghost layer.
        algo = espressomd.reaction_methods.ExclusionRadius(
            exclusion_range=0.5, search_algorithm="order_n")

        # inserted particle == particle 0; particle 1 lies within the exclusion
        # range and must be detected.
        self.assertTrue(algo.check_exclusion_range(pid=0, ptype=0))
        self.assertTrue(algo.check_exclusion_range(pid=0, ptype=-2))
        self.assertTrue(algo.check_exclusion_range(pid=0))

        # symmetric case: inserted particle == particle 1
        self.assertTrue(algo.check_exclusion_range(pid=1, ptype=0))
        self.assertTrue(algo.check_exclusion_range(pid=1, ptype=-2))
        self.assertTrue(algo.check_exclusion_range(pid=1))

        # negative control: a small exclusion range below the mi-distance must
        # NOT flag an overlap (the order_n branch must not over-detect either)
        algo = espressomd.reaction_methods.ExclusionRadius(
            exclusion_range=0.05, search_algorithm="order_n")
        self.assertFalse(algo.check_exclusion_range(pid=0, ptype=0))
        self.assertFalse(algo.check_exclusion_range(pid=0, ptype=-2))
        self.assertFalse(algo.check_exclusion_range(pid=0))


if __name__ == "__main__":
    ut.main()
