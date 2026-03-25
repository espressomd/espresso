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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.polymer
import espressomd.interactions
import numpy as np
import contextlib
import io


@utx.skipIfMissingFeatures(["VIRTUAL_SITES_CENTER_OF_MASS"])
class VirtualSitesCOM(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    FENE_PARAMS = {'k': 7, 'r_0': 1, 'd_r_max': 2}
    fene = espressomd.interactions.FeneBond(**FENE_PARAMS)
    POLYMER_PARAMS = {'n_polymers': 1, 'bond_length': 1,
                      'seed': 42, 'min_distance': 0.9}

    np.random.seed(42)

    def build_polymer(self, n_monomers, polymer_params,
                      fene, monomer_type=0, mol_id=0):
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

    def set_molecules_and_vs(
            self, molecule_ids, n_monomers, monomer_types, vs_type=5):
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
        for molecule_id, n_monomers, monomer_type in zip(
                molecule_ids, n_monomers, monomer_types):
            # Build polymer chain
            self.build_polymer(n_monomers, self.POLYMER_PARAMS,
                               self.fene, monomer_type, molecule_id)
            # Add virtual particle at the origin
            vs = self.system.part.add(pos=[0, 0, 0], type=vs_type)
            vs.vs_com_relate_to(molecule_id)
            mid_for_vs[molecule_id] = vs.id

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
        self.system.bonded_inter.clear()

    def test_vs_initialization(self):
        """
        Test initialization of virtual sites com
        """

        p1 = self.system.part.add(
            pos=[0, 0, 0], virtual=True, type=1, id=1, mol_id=10)
        vs1 = self.system.part.add(pos=[1, 1, 1], virtual=True, type=1, id=2)
        vs1.vs_com_relate_to(p1)
        self.assertEqual(vs1.mol_id, p1.mol_id)

        mol_id_p2 = 20
        p2 = self.system.part.add(
            pos=[2, 2, 2], virtual=True, type=1, id=3, mol_id=mol_id_p2)
        vs2 = self.system.part.add(pos=[3, 3, 3], virtual=True, type=1, id=4)
        vs2.vs_com_relate_to(mol_id_p2)
        self.assertEqual(vs2.mol_id, p2.mol_id)

    def test_vs_position_mass(self):
        """
        Test update of the vs positions and masses
        """

        molecule_ids = [1, 2]
        n_monomers = [20, 50]
        monomer_types = [0, 1]

        mid_for_vs = self.set_molecules_and_vs(
            molecule_ids, n_monomers, monomer_types)

        self.system.integrator.set_steepest_descent(
            f_max=10., gamma=1., max_displacement=0.2)
        self.system.integrator.run(100)

        # Check position of virtual sites after a steepest descent integration
        for mol_id, vs_id, monomer_type in zip(
                mid_for_vs.keys(), mid_for_vs.values(), monomer_types):
            # test vs position
            vs_pos = self.system.part.by_id(vs_id).pos
            expected_vs_pos = self.system.analysis.center_of_mass(
                p_type=monomer_type)
            for pair in zip(expected_vs_pos, vs_pos):
                self.assertAlmostEqual(pair[0], pair[1])
            # test vs mass
            vs_mass = self.system.part.by_id(vs_id).mass
            expected_vs_mass = 0.
            for p in self.system.part.select(mol_id=mol_id):
                if not p.is_virtual():
                    expected_vs_mass += p.mass
            self.assertEqual(vs_mass, expected_vs_mass)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_particle_forces(self):
        """
        Test force on molecule particles when the vs undergoes external force.
        """

        molecule_id = [1]
        # larger polymer to have a complete distribution over the ranks
        n_monomers = [200]
        monomer_types = [0]
        applied_force = np.array([100., 0., 0.])

        mid_for_vs = self.set_molecules_and_vs(
            molecule_id, n_monomers, monomer_types)

        self.system.integrator.set_steepest_descent(
            f_max=10., gamma=1., max_displacement=0.2)
        self.system.integrator.run(100)

        vs_part = self.system.part.by_id(mid_for_vs[molecule_id[0]])
        vs_part.ext_force = applied_force
        ref_force = applied_force / n_monomers[0]

        self.system.integrator.run(1)

        for p in self.system.part.select(mol_id=molecule_id[0]):
            if not p.is_virtual():
                np.testing.assert_allclose(np.copy(p.f), ref_force, atol=1e-7)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_umbrella_sampling(self):
        """
        Check umbrella sampling. Introduce a harmonic bias to keep the polymer
        close to a reference point that is constrained to move on a line.
        """

        molecule_id = [1]
        n_monomers = [200]
        monomer_types = [0]
        harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0., k=1e4)
        self.system.bonded_inter.add(harmonic_bond)

        mid_for_vs = self.set_molecules_and_vs(
            molecule_id, n_monomers, monomer_types)

        self.system.integrator.set_steepest_descent(
            f_max=10., gamma=1., max_displacement=0.2)
        self.system.integrator.run(100)
        self.system.integrator.set_vv()
        self.system.thermostat.set_langevin(kT=1., gamma=1., seed=1)

        # tether VS COM to a particle constrained to move on a reaction path
        vs_part = self.system.part.by_id(mid_for_vs[molecule_id[0]])
        p_tether = self.system.part.add(pos=vs_part.pos + 1e-6,
                                        type=5, fix=[True, True, False])
        vs_part.add_bond((harmonic_bond, p_tether))

        # track position of the center of mass over time
        xdata = np.zeros([3000, 3])
        for i in range(xdata.shape[0]):
            self.system.integrator.run(4)
            xdata[i] = self.system.analysis.center_of_mass(p_type=0)

        # the standard deviation of a constrained coordinate is directly
        # related to the strength of the harmonic potential
        self.assertAlmostEqual(np.std(xdata[:, 0]), 10e-3, delta=3e-3)
        self.assertAlmostEqual(np.std(xdata[:, 1]), 10e-3, delta=3e-3)
        # the standard deviation of the other coordinates is not bounded
        self.assertGreater(np.std(xdata[:, 2]), 0.1)

    def test_vs_exceptions(self):
        """
        Test exceptions related to virtual sites com
        """
        p = self.system.part.add(pos=0.9 * self.system.box_l, type=1, mol_id=8)
        vs1 = self.system.part.add(pos=[0, 0, 0], type=1)
        vs2 = self.system.part.add(pos=[1, 1, 1], type=1)

        vs1.vs_com_relate_to(p)
        # relate to empty
        with self.assertRaisesRegex(TypeError, "missing 1 required positional argument"):
            vs1.vs_com_relate_to()
        # relating to anything else other than a particle or id is not allowed
        with self.assertRaisesRegex(ValueError, "Argument of 'vs_com_relate_to' has to be of type ParticleHandle or int"):
            vs1.vs_com_relate_to('0')
        with self.assertRaisesRegex(ValueError, "Invalid molecule id: -2"):
            vs1.vs_com_relate_to(-2)
        # relating to the same molecule twice is not allowed
        with self.assertRaisesRegex(RuntimeError, "Molecule id: 8 is already tracked by virtual site with particle id: 1"):
            vs2.vs_com_relate_to(p)
        # state remains unchanged
        self.assertEqual(vs1.mol_id, p.mol_id)
        self.assertTrue(vs1.is_virtual())
        self.assertFalse(vs2.is_virtual())

        with contextlib.redirect_stderr(io.StringIO()):
            if espressomd.has_features("EXTERNAL_FORCES"):
                vs1.fix = [True, False, False]
                with self.assertRaisesRegex(Exception, "VS COM particles cannot be fixed in space"):
                    self.system.integrator.run(0)
                vs1.fix = 3 * [False]
            if espressomd.has_features("BOND_CONSTRAINT"):
                rigid_bond = espressomd.interactions.RigidBond(
                    r=0.1, ptol=1E-2, vtol=1e-2)
                self.system.bonded_inter.add(rigid_bond)
                vs1.add_bond((rigid_bond, p))
                with self.assertRaisesRegex(Exception, "VS COM particles cannot use rigid bonds"):
                    self.system.integrator.run(0)
                vs1.delete_bond((rigid_bond, p))
        # system can still be integrated
        self.system.integrator.run(0)

        # relating to own molecule is allowed
        vs3 = self.system.part.add(pos=[0, 0, 0], mol_id=p.mol_id + 1)
        vs3.vs_com_relate_to(vs3.mol_id)
        self.assertEqual(vs3.mol_id, p.mol_id + 1)


if __name__ == "__main__":
    ut.main()
