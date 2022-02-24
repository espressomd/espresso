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

import espressomd
import espressomd.reaction_ensemble

import unittest as ut


class ReactionMethods(ut.TestCase):

    """Test the reaction methods interface."""

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.part.clear()

    def check_interface(self, method, kT, exclusion_radius, gamma):
        def check_reaction_parameters(reactions, parameters):
            for reaction, params in zip(reactions, parameters):
                for key in reaction.required_keys():
                    self.assertEqual(getattr(reaction, key), params[key])

        reaction_forward = {
            'gamma': gamma,
            'reactant_types': [5],
            'reactant_coefficients': [1],
            'product_types': [2, 3],
            'product_coefficients': [1, 1],
            'default_charges': {5: 0, 2: 0, 3: 0},
        }
        reaction_backward = {
            'gamma': 1. / gamma,
            'reactant_types': reaction_forward['product_types'],
            'reactant_coefficients': reaction_forward['product_coefficients'],
            'product_types': reaction_forward['reactant_types'],
            'product_coefficients': reaction_forward['reactant_coefficients'],
            'default_charges': reaction_forward['default_charges'],
        }

        if isinstance(method, espressomd.reaction_ensemble.ConstantpHEnsemble):
            method.add_reaction(gamma=reaction_forward['gamma'],
                                reactant_types=reaction_forward['reactant_types'],
                                product_types=reaction_forward['product_types'],
                                default_charges=reaction_forward['default_charges'])
        else:
            method.add_reaction(**reaction_forward)
        reaction_parameters = (reaction_forward, reaction_backward)

        # check getters and setters
        self.assertAlmostEqual(method.kT, kT, delta=1e-10)
        self.assertAlmostEqual(
            method.exclusion_radius,
            exclusion_radius,
            delta=1e-10)
        self.assertAlmostEqual(
            method.get_volume(),
            self.system.volume(),
            delta=1e-10)
        method.set_volume(volume=1.)
        self.assertAlmostEqual(method.get_volume(), 1., delta=1e-10)
        self.assertEqual(method.get_non_interacting_type(), 100)
        method.set_non_interacting_type(type=9)
        self.assertEqual(method.get_non_interacting_type(), 9)
        if isinstance(method, espressomd.reaction_ensemble.ConstantpHEnsemble):
            self.assertAlmostEqual(method.constant_pH, 10., delta=1e-10)
            method.constant_pH = 8.
            self.assertAlmostEqual(method.constant_pH, 8., delta=1e-10)

        # check constraints
        method.set_wall_constraints_in_z_direction(
            slab_start_z=0.1, slab_end_z=0.9)
        offsets = method.get_wall_constraints_in_z_direction()
        self.assertAlmostEqual(offsets[0], 0.1, delta=1e-10)
        self.assertAlmostEqual(offsets[1], 0.9, delta=1e-10)
        method.remove_constraint()

        # check status
        status = method.get_status()
        self.assertEqual(status['kT'], kT)
        self.assertEqual(status['exclusion_radius'], exclusion_radius)
        self.assertEqual(len(status['reactions']), 2)
        for reaction_flat, params in zip(
                status['reactions'], reaction_parameters):
            for key in reaction_flat:
                if key == 'gamma':
                    self.assertAlmostEqual(
                        reaction_flat[key], params[key], delta=1e-10)
                else:
                    self.assertEqual(reaction_flat[key], params[key])

        # check reactions
        reactions = method.reactions
        self.assertEqual(len(reactions), 2)
        check_reaction_parameters(method.reactions, reaction_parameters)

        # check reactions after parameter change
        new_gamma = 634.
        reaction_forward['gamma'] = new_gamma
        reaction_backward['gamma'] = 1. / new_gamma
        method.change_reaction_constant(reaction_id=0, gamma=new_gamma)
        check_reaction_parameters(method.reactions, reaction_parameters)
        status = method.get_status()
        self.assertAlmostEqual(
            status['reactions'][0]['gamma'],
            reaction_forward['gamma'],
            delta=1e-10)
        self.assertAlmostEqual(
            status['reactions'][1]['gamma'],
            reaction_backward['gamma'],
            delta=1e-10)

        # check particle deletion
        p1, _, p3 = self.system.part.add(
            pos=3 * [(0., 0., 0.)], type=[5, 2, 3])
        if isinstance(method, espressomd.reaction_ensemble.WidomInsertion):
            potential_energy = method.calculate_particle_insertion_potential_energy(
                reaction_id=0)
            self.assertEqual(potential_energy, 0.)
        method.delete_particle(p_id=p3.id)
        self.assertEqual(len(self.system.part), 2)
        method.delete_particle(p_id=p1.id)
        self.assertEqual(len(self.system.part), 1)
        self.system.part.clear()

        # check reaction deletion
        method.delete_reaction(reaction_id=0)
        self.assertEqual(len(method.reactions), 0)

    def test_interface(self):
        # reaction ensemble
        method = espressomd.reaction_ensemble.ReactionEnsemble(
            kT=1.5, exclusion_radius=0.8, seed=12)
        self.check_interface(method, kT=1.5, exclusion_radius=0.8, gamma=1.2)

        # constant pH ensemble
        method = espressomd.reaction_ensemble.ConstantpHEnsemble(
            kT=1.5, exclusion_radius=0.8, seed=12, constant_pH=10)
        self.check_interface(method, kT=1.5, exclusion_radius=0.8, gamma=1.2)

        # Widom insertion
        method = espressomd.reaction_ensemble.WidomInsertion(kT=1.6, seed=12)
        self.check_interface(method, kT=1.6, exclusion_radius=0., gamma=1.)

    def test_exceptions(self):
        single_reaction_params = {
            'gamma': 1.,
            'reactant_types': [4],
            'reactant_coefficients': [1],
            'product_types': [2, 3],
            'product_coefficients': [1, 4],
        }
        reaction_params = {
            'default_charges': {2: 0, 3: 0, 4: 0},
            **single_reaction_params
        }
        widom = espressomd.reaction_ensemble.WidomInsertion(kT=1., seed=12)
        method = espressomd.reaction_ensemble.ReactionEnsemble(
            kT=1.5, exclusion_radius=0.8, seed=12)
        method.add_reaction(**reaction_params)
        widom.add_reaction(**reaction_params)

        # check invalid reactions
        err_msg = 'number of types and coefficients have to match'
        with self.assertRaisesRegex(ValueError, f'reactants: {err_msg}'):
            method.add_reaction(**{**reaction_params, 'reactant_types': []})
        with self.assertRaisesRegex(ValueError, f'products: {err_msg}'):
            method.add_reaction(**{**reaction_params, 'product_types': []})

        # check charge conservation
        err_msg = 'Reaction system is not charge neutral'
        with self.assertRaisesRegex(ValueError, err_msg):
            method.add_reaction(default_charges={2: 8, 3: 0, 4: -50},
                                **single_reaction_params)
        with self.assertRaisesRegex(ValueError, err_msg):
            method.add_reaction(default_charges={2: 1, 3: 0, 4: 1 + 1e-10},
                                **single_reaction_params)

        # check invalid reaction id exceptions
        # (note: reactions id = 2 * reactions index)
        self.assertEqual(len(method.reactions), 2)
        for i in [-2, -1, 1, 2, 3]:
            with self.assertRaisesRegex(IndexError, 'This reaction is not present'):
                method.delete_reaction(reaction_id=i)
            with self.assertRaisesRegex(IndexError, 'This reaction is not present'):
                method.get_acceptance_rate_reaction(reaction_id=2 * i)

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
            method.delete_particle(p_id=0)
        with self.assertRaisesRegex(RuntimeError, "Trying to remove some non-existing particles from the system via the inverse Widom scheme"):
            widom.calculate_particle_insertion_potential_energy(reaction_id=0)

        # check other exceptions
        with self.assertRaisesRegex(ValueError, "Invalid value for 'volume'"):
            method.set_volume(volume=-10.)
        with self.assertRaisesRegex(RuntimeError, r"unknown method 'unknown\(\)'"):
            method.call_method('unknown', x=1)
        err_msg = r"Only the following keys can be given as keyword arguments: \[.+\], got \[.+\] \(unknown \['x'\]\)"
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_ensemble.SingleReaction(
                x=1, **single_reaction_params)
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_ensemble.ReactionEnsemble(
                kT=1., exclusion_radius=1., seed=12, x=1)
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_ensemble.ConstantpHEnsemble(
                kT=1., exclusion_radius=1., seed=12, x=1, constant_pH=2)
        with self.assertRaisesRegex(ValueError, err_msg):
            espressomd.reaction_ensemble.WidomInsertion(
                kT=1., seed=12, x=1)
        with self.assertRaisesRegex(ValueError, "Invalid value for 'kT'"):
            espressomd.reaction_ensemble.ReactionEnsemble(
                kT=-1., exclusion_radius=1., seed=12)
        with self.assertRaisesRegex(ValueError, "Invalid value for 'exclusion_radius'"):
            espressomd.reaction_ensemble.ReactionEnsemble(
                kT=1., exclusion_radius=-1., seed=12)


if __name__ == "__main__":
    ut.main()
