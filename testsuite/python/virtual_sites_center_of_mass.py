#
# Copyright (C) 2013-2025 The ESPResSo project
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
import espressomd.polymer
import numpy as np


@utx.skipIfMissingFeatures(["VIRTUAL_SITES_CENTER_OF_MASS"])
class VirtualSitesCOM(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    FENE_PARAMS = {'k': 7, 'r_0': 1, 'd_r_max': 2}
    fene = espressomd.interactions.FeneBond(**FENE_PARAMS)
    POLYMER_PARAMS = {'n_polymers': 1, 'bond_length': 1, 'seed': 42, 'min_distance': 0.9}
    
    np.random.seed(42)

    def build_polymer(self, n_monomers, polymer_params, fene, monomer_type=0, mol_id=0):
        """
        Build a polymer chain with the specified number of monomers, bond type
        and molecule id.
        """
        positions = espressomd.polymer.linear_polymer_positions(
            beads_per_chain=n_monomers, **polymer_params)
        p_previous = None
        for pos in positions[0]:
            p = self.system.part.add(pos=pos, mol_id=mol_id, type=monomer_type)
            if p_previous is not None:
                p.add_bond((fene, p_previous))
            p_previous = p

    def set_molecules_and_vs(self, molecule_ids, n_monomers, monomer_types, vs_type=5, id_shift=10):
        """
        Set virtual sites and the corresponding polymer molecules.

        This function creates polymer molecules with the specified number of monomers
        and types, and the corresponding virtual site with each molecule. The virtual site
        is added at the origin and is assigned a unique ID.

        Parameters:
        - molecule_ids (list of int): List of molecule IDs of the polymer chains.
        - n_monomers (list of int): List of the number of monomers in each polymer chain.
        - monomer_types (list of int): List of the types of monomers in each polymer chain.

        Returns:
        - dict: A dictionary mapping molecule IDs to the IDs of their corresponding virtual sites.
        """
        mid_for_vs = {}
        for molecule_id_,n_monomers_,monomer_type_ in zip(molecule_ids, n_monomers, monomer_types):
            # Build polymer chain
            self.build_polymer(n_monomers_, self.POLYMER_PARAMS, self.fene, monomer_type_, molecule_id_)
            # Add virtual particle at the origin
            vs = self.system.part.add(pos=[0, 0, 0], virtual=True, type=vs_type, mol_id=molecule_id_+id_shift)
            vs.vs_com_auto_relate_to(molecule_id_)
            mid_for_vs[molecule_id_] = vs.id

        return mid_for_vs

    def setUp(self):
        self.system.box_l = [30.0, 30.0, 30.0]
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4
        self.system.bonded_inter.add(self.fene)

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()


    def test_vs_position_mass(self):
        """
        Test update of the vs positions and masses
        """

        molecule_ids = [1, 2]
        n_monomers = [20, 50]
        monomer_types = [0, 1]

        mid_for_vs = self.set_molecules_and_vs(molecule_ids, n_monomers, monomer_types)

        self.system.integrator.set_steepest_descent(f_max=10, gamma=50.0, max_displacement=0.2)
        self.system.integrator.run(1)

        # Check position of virtual sites after a steepest descent inegration
        for mol_id_,vs_id_,monomer_type_ in zip(mid_for_vs.keys(), mid_for_vs.values(), monomer_types):
            # test vs position
            vs_pos = self.system.part.by_id(vs_id_).pos
            expected_vs_pos = self.system.analysis.center_of_mass(p_type=monomer_type_)
            for pair in zip(expected_vs_pos, vs_pos):
                self.assertAlmostEqual(pair[0], pair[1])
            # test vs mass
            vs_mass = self.system.part.by_id(vs_id_).mass
            expected_vs_mass = 0
            for part in self.system.part.select(mol_id=mol_id_):
                expected_vs_mass += part.mass
            self.assertEqual(expected_vs_mass, vs_mass)


    # def test_particle_forces(self):
    #     """
    #     Test force on molecule particles when the vs undergoes given force 
    #     """

    #     molecule_id = [1]
    #     n_monomers = [50]
    #     monomer_types = [0]
    #     applied_force = np.array([100, 0, 0], dtype=float)

    #     mid_for_vs = self.set_molecules_and_vs(molecule_id, n_monomers, monomer_types)

    #     self.system.integrator.set_steepest_descent(f_max=10, gamma=50.0, max_displacement=0.2)
    #     self.system.integrator.run(1000)

    #     vs_part = self.system.part.by_id(mid_for_vs[molecule_id[0]])
    #     vs_part.ext_force = applied_force
    #     expected_force = applied_force/n_monomers[0]
        
    #     self.system.integrator.run(1)

    #     for part in self.system.part.select(mol_id=molecule_id[0]):
    #         for pair in zip(expected_force, part.f):
    #             self.assertAlmostEqual(pair[0], pair[1])


    # def test_vs_exceptions(self):
    #     """
    #     Test exceptions related to virtual sites com
    #     """
    #     vs1 = self.system.part.add(pos=[0, 0, 0], virtual=True, type=1, id=1)
    #     vs2 = self.system.part.add(pos=[1, 1, 1], virtual=True, type=1, id=2)
    #     # relate to empty
    #     with self.assertRaisesRegex(TypeError, "Argument of 'vs_com_auto_relate_to' has to be of type int list of integers"):
    #         vs1.vs_com_auto_relate_to()
    #     # relating to anything else other than a particle or id is not allowed
    #     with self.assertRaisesRegex(ValueError, "Argument of 'vs_com_auto_relate_to' has to be of type int list of integers"):
    #         vs1.vs_com_auto_relate_to('0')
    #     with self.assertRaisesRegex(ValueError, "Invalid particle id: -2"):
    #         vs1.vs_com_auto_relate_to(-2)
    #     # relating to itself is not allowed
    #     with self.assertRaisesRegex(ValueError, "A virtual site cannot relate to itself"):
    #         vs1.vs_com_auto_relate_to(vs1)
    #     # relating to a non-existing particle id is not allowed
    #     with self.assertRaisesRegex(ValueError, "No real particle with id 3 for virtual site with id 1"):
    #         vs1.vs_com_auto_relate_to(3)










    #     # dangling virtual sites are not allowed
    #     with self.assertRaisesRegex(Exception, "Particle with id 4 is a dangling virtual site"):
    #         p4.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_RELATIVE
    #         self.assertEqual(p4.vs_relative[0], -1)
    #         system.integrator.run(0, recalc_forces=True)
    #     p4.remove()
    #     # relating to a deleted particle is not allowed
    #     with self.assertRaisesRegex(Exception, "No real particle with id 3 for virtual site with id 2"):
    #         p2.vs_auto_relate_to(p3)
    #         p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_RELATIVE
    #         p3.remove()
    #         system.integrator.run(0, recalc_forces=True)
    #     if system.cell_system.get_state()["n_nodes"] > 1:
    #         with self.assertRaisesRegex(Exception, r"The distance between virtual and non-virtual particle \([0-9\.]+\) is larger than the minimum global cutoff"):
    #             p2.vs_auto_relate_to(p1)
    #         # If overridden this check should not raise an exception
    #         p2.vs_auto_relate_to(p1, override_cutoff_check=True)

    # def test_exceptions(self):
    #     """
    #     """


    #     self.assertRaises(vs.vs_com_auto_relate_to([]))


if __name__ == "__main__":
    ut.main(verbosity=2)

# virtual site is a particle with its properties -> vs have also mol_id=0 by default 
# Naming: vs.vs_com_auto_relate_to() ok?