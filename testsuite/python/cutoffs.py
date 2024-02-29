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
from espressomd.interactions import FeneBond
import numpy as np
import unittest as ut
import unittest_decorators as utx


class CutOff(ut.TestCase):

    """
    Test interaction cutoffs.
    Tests must be executed in the order given below such 
    that the leftovers from test 2 do not break test 1.
    """
    system = espressomd.System(box_l=3 * [50])

    def tearDown(self) -> None:
        self.system.non_bonded_inter.reset()
        self.system.bonded_inter.clear()
        self.system.part.clear()

    def test_1max_cut(self):
        system = self.system
        system.cell_system.skin = 1

        # Initial state. Skin does not influence cutoffs as long as there are
        # no interactions
        self.assertEqual(system.cell_system.max_cut_nonbonded, -1)
        self.assertEqual(system.cell_system.max_cut_bonded, -1)
        self.assertEqual(system.cell_system.interaction_range, -1)

        # Bonded interaction
        fene = FeneBond(r_0=1, d_r_max=2, k=1)
        system.bonded_inter.add(fene)
        self.assertEqual(system.cell_system.max_cut_bonded, 3)
        n_nodes = np.product(system.cell_system.node_grid)
        if n_nodes == 1:
            # Bonds don't influence interaction range
            self.assertEqual(system.cell_system.interaction_range, -1)
        else:
            self.assertEqual(system.cell_system.interaction_range,
                             system.cell_system.max_cut_bonded +
                             system.cell_system.skin)

        system.bonded_inter.remove(fene._bond_id)
        self.assertEqual(system.cell_system.max_cut_bonded, -1)
        self.assertEqual(system.cell_system.interaction_range, -1)

        if espressomd.has_features("LENNARD_JONES"):
            system.non_bonded_inter[0, 0].lennard_jones.set_params(
                sigma=1, epsilon=1, cutoff=2.5, shift="auto")
            self.assertEqual(
                system.cell_system.max_cut_nonbonded, 2.5)
            self.assertEqual(system.cell_system.interaction_range,
                             system.cell_system.max_cut_nonbonded +
                             system.cell_system.skin)

            system.non_bonded_inter[0, 0].lennard_jones.deactivate()
            self.assertEqual(system.cell_system.max_cut_nonbonded, -1)
            self.assertEqual(system.cell_system.interaction_range, -1)

    @utx.skipIfMissingFeatures(["LENNARD_JONES",
                               "VIRTUAL_SITES_RELATIVE", "COLLISION_DETECTION", "DP3M"])
    def test_2cutoff_by_type(self):
        import espressomd.virtual_sites

        sys = self.system
        min_global_cut = 0.1
        sys.min_global_cut = min_global_cut
        sys.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=0.1, epsilon=1, cutoff=1, shift="auto")
        sys.non_bonded_inter[0, 1].lennard_jones.set_params(
            sigma=0.1, epsilon=1, cutoff=4, shift="auto")
        sys.non_bonded_inter[0, 2].lennard_jones.set_params(
            sigma=0.1, epsilon=1, cutoff=3, shift="auto")

        self.assertEqual(sys.max_cut_nonbonded, 4)
        self.assertEqual(sys.cutoff_by_types([0]), 1)
        self.assertEqual(sys.cutoff_by_types([0, 1, 2]), 4)
        self.assertEqual(sys.cutoff_by_types([0, 2, 100]), 3)
        self.assertEqual(sys.cutoff_by_types([1, 2]), min_global_cut)
        self.assertEqual(sys.cutoff_by_types(
            [9, 10, 11, 23, 1456]), min_global_cut)

        sys.non_bonded_inter[0, 2].lennard_jones.set_params(
            sigma=0.1, epsilon=1, cutoff=2, shift="auto")
        self.assertEqual(sys.cutoff_by_types([0, 2]), 2)

        sys.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
        bond = espressomd.interactions.HarmonicBond(k=1000, r_0=0)
        sys.bonded_inter.add(bond)
        sys.collision_detection.set_params(
            mode="glue_to_surface",
            distance=5.5,
            distance_glued_particle_to_vs=0.02,
            bond_centers=bond,
            bond_vs=bond,
            part_type_vs=1,
            part_type_to_attach_vs_to=2,
            part_type_to_be_glued=7,
            part_type_after_glueing=8)
        self.assertEqual(sys.cutoff_by_types([2, 7]), 5.5)
        self.assertEqual(sys.cutoff_by_types([1, 7]), min_global_cut)
        self.assertEqual(sys.cutoff_by_types([7, 8]), min_global_cut)

        sys.collision_detection.set_params(mode="bind_centers", distance=6.6,
                                           bond_centers=bond)
        self.assertEqual(sys.cutoff_by_types([1, 2]), 6.6)
        self.assertEqual(sys.cutoff_by_types([1, 7]), 6.6)
        self.assertEqual(sys.cutoff_by_types([7, 8]), 6.6)
        import espressomd.magnetostatics as magnetostatics
        p3m = magnetostatics.DipolarP3M(
            prefactor=1,
            mesh=32,
            alpha=0.123,
            accuracy=1E-4,
            r_cut=7.7,
            cao=3,
            tune=False)
        sys.part.add(pos=50 * np.random.random((10, 3)),
                     dip=np.random.random((10, 3)))
        sys.magnetostatics.solver = p3m
        self.assertEqual(sys.cutoff_by_types([1, 2]), 7.7)


if __name__ == '__main__':
    ut.main()
