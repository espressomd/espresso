# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import numpy as np
from espressomd import System, polymer, interactions


class DiamondPolymer(ut.TestCase):
    """
    Test the functionality of espressomd.polymer.setup_diamond_polymer()
    in terms of
    * properties of particles created
    * connections via bonds
    * the geometry of the polymer network
    """

    system = System(box_l=3 * [16])
    diamond_params = {'MPC': 15, 
                      'dist_cM': 3, 
                      'val_cM': -1.3, 
                      'val_nodes': 0.6, 
                      'start_id': 3, 
                      'no_bonds': False, 
                      'type_nodes': 2,
                      'type_nM': 5, 
                      'type_cM': 7}
    bond_length = system.box_l[0] * \
        (0.25 * np.sqrt(3)) / (diamond_params['MPC'] + 1)

    def setUp(self):
        bond = interactions.HarmonicBond(k=1.5, r_0=self.bond_length, r_cut=3) 
        self.system.bonded_inter.add(bond)
        polymer.setup_diamond_polymer(system=self.system,
                                      bond=bond,
                                      **self.diamond_params)
        self.system.time_step = 0.1
        self.node_parts = self.system.part.select(
            type=self.diamond_params['type_nodes'])
        self.charged_mono = self.system.part.select(
            type=self.diamond_params['type_cM'])
        self.noncharged_mono = self.system.part.select(
            type=self.diamond_params['type_nM'])

    def tearDown(self):
        self.system.part.clear()

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_particle_properties(self):
        """
        checks if the particles created have the right type and charge

        """    
        # number of particles created
        number_non_node = 16 * self.diamond_params['MPC']
        self.assertEqual(len(self.system.part), number_non_node + 8)
        # number of each type
        self.assertEqual(len(self.node_parts), 8)
        self.assertEqual(len(self.charged_mono), 
                         np.rint(number_non_node / self.diamond_params['dist_cM']))
        self.assertEqual(len(self.noncharged_mono),
                         np.rint(number_non_node * (1 - 1 / self.diamond_params['dist_cM'])))
        # charge
        np.testing.assert_allclose(self.node_parts.q, 
                                   len(self.node_parts) * [self.diamond_params['val_nodes']])
        np.testing.assert_allclose(self.noncharged_mono.q,
                                   len(self.noncharged_mono) * [0])                                                                 
        np.testing.assert_allclose(self.charged_mono.q, 
                                   len(self.charged_mono) * [self.diamond_params['val_cM']])
        # particle id
        self.assertGreaterEqual(min(self.system.part[:].id), 
                                self.diamond_params['start_id'])

    def test_bonds(self):
        """
        test that the right number of bonds is formed on each particle
        """  
        # 4 bonds on nodes
        for part in self.node_parts:
            self.assertEqual(len(part.bonds), 4)
        # 1 or 0 bonds on chain monomers
        n_bonds = [len(bonds) for bonds in 
                   self.noncharged_mono.bonds + self.charged_mono.bonds]
        self.assertLessEqual(np.max(n_bonds), 1)
        # total number of bonds
        number_non_node = 16 * self.diamond_params['MPC']
        self.assertEqual(np.sum(n_bonds), number_non_node - 16)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_connected(self):
        """
        test that all particles in the polymer are connected by pushing one particle
        """

        self.system.part[self.diamond_params['start_id']].ext_force = 3 * [2]
        self.system.integrator.run(200)
        vels = np.linalg.norm(self.system.part[:].v, axis=1)
        self.assertGreater(np.min(vels), 1e-6)

    def test_geometry(self):
        """
        check if the distance between all monomers is correct,
        only nearest neighbouring nodes are connected
        and that the nodes have the right position

        """
        # Energy calculation checks distance indirectly through
        # position of minimum of HarmonicBond
        # With formula for self.bond_length this also ensures
        # that only nearest neighbours can be reached
        E = self.system.analysis.energy()['total']
        self.assertAlmostEqual(E, 0., delta=1e-13)

        node_pos_scaled = np.array(self.node_parts.pos) / self.system.box_l[0]
        node_pos_shouldbe = 0.25 * np.array([[0, 0, 0], [1, 1, 1],
                                             [2, 2, 0], [0, 2, 2],
                                             [2, 0, 2], [3, 3, 1],
                                             [1, 3, 3], [3, 1, 3]])
        for pos1, pos2 in zip(node_pos_scaled, node_pos_shouldbe):
            np.testing.assert_allclose(pos1, pos2)


if __name__ == "__main__":
    ut.main()
