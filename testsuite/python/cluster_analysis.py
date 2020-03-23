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
import espressomd
from espressomd.utils import handle_errors
import numpy as np
from espressomd.interactions import FeneBond
from espressomd.pair_criteria import DistanceCriterion, BondCriterion
from espressomd.cluster_analysis import ClusterStructure


class ClusterAnalysis(ut.TestCase):

    """Tests the cluster analysis"""

    es = espressomd.System(box_l=(1, 1, 1))

    f = FeneBond(k=1, d_r_max=0.05)
    es.bonded_inter.add(f)

    # 1st cluster
    es.part.add(id=0, pos=(0, 0, 0))
    es.part.add(id=1, pos=(0.91, 0, 0), bonds=((0, 0),))
    es.part.add(id=2, pos=(0, 0.2, 0))
    es.part.add(id=3, pos=(0, 0.1, 0))

    # 2nd cluster
    es.part.add(id=4, pos=(0.5, 0.5, 0.5))
    es.part.add(id=5, pos=(0.55, 0.5, 0.5))

    cs = ClusterStructure()
    np.random.seed(1)

    # Setup check
    handle_errors("")

    def test_00_fails_without_criterion_set(self):
        with self.assertRaises(Exception):
            self.cs.run_for_all_pairs()

    def test_set_criterion(self):
        # Test setters/getters for criteria
        dc = DistanceCriterion(cut_off=0.11)
        self.cs.set_params(pair_criterion=dc)
        # Do we get back the right criterion
        dc_ret = self.cs.get_params()["pair_criterion"]
        # Note: This work around the fact that the script interface does not
        # yet assign the correct derived class when returning an object
        self.assertEqual(dc_ret.name(), "PairCriteria::DistanceCriterion")
        self.assertAlmostEqual(dc_ret.get_params()["cut_off"], 0.11, places=7)

        # Is the cluster structure empty before being used

    def test_analysis_for_all_pairs(self):
        # Run cluster analysis
        self.cs.set_params(pair_criterion=DistanceCriterion(cut_off=0.12))
        self.cs.run_for_all_pairs()

        # Number of clusters
        self.assertEqual(len(self.cs.clusters), 2)
        cids = self.cs.cluster_ids()

        # Sizes of individual clusters
        l2 = self.cs.clusters[cids[1]].size()
        l1 = len(self.cs.clusters[cids[0]].particle_ids())

        # Clusters should contain 2 and 4 particles
        self.assertEqual(min(l1, l2), 2)
        self.assertEqual(max(l1, l2), 4)

        # Verify particle ids
        smaller_cluster = None
        bigger_cluster = None
        if l1 < l2:
            smaller_cluster = self.cs.clusters[cids[0]]
            bigger_cluster = self.cs.clusters[cids[1]]
        else:
            smaller_cluster = self.cs.clusters[cids[1]]
            bigger_cluster = self.cs.clusters[cids[0]]

        self.assertEqual(bigger_cluster.particle_ids(), [0, 1, 2, 3])
        self.assertEqual(smaller_cluster.particle_ids(), [4, 5])

        # Test obtaining a ParticleSlice for a cluster
        pids = bigger_cluster.particle_ids()
        particles = bigger_cluster.particles()
        # Do the number of entries match
        self.assertEqual(len(pids), len(particles.id_selection))

        # Compare ids of particles in the slice
        self.assertEqual(all(particles.id_selection), all(pids))

        # Test iteration over clusters
        visited_sizes = []
        for c in self.cs.clusters:
            visited_sizes.append(c[1].size())
        visited_sizes = sorted(visited_sizes)
        self.assertEqual(visited_sizes, [2, 4])

    def test_zz_single_cluster_analysis(self):
        self.es.part.clear()
        # Place particles on a line (crossing periodic boundaries)
        for x in np.arange(-0.2, 0.21, 0.01):
            self.es.part.add(pos=(x, 1.1 * x, 1.2 * x))
        self.cs.pair_criterion = DistanceCriterion(cut_off=0.13)
        self.cs.run_for_all_pairs()
        self.assertEqual(len(self.cs.clusters), 1)

        for c in self.cs.clusters:
            # Discard cluster id
            c = c[1]

            # Center of mass should be at origin
            self.assertLess(np.linalg.norm(c.center_of_mass()), 1E-8)

            # Longest distance
            self.assertAlmostEqual(
                c.longest_distance(),
                self.es.distance(self.es.part[0],
                                 self.es.part[len(self.es.part) - 1]),
                delta=1E-8)

            # Radius of gyration
            rg = 0.
            com_particle = self.es.part[len(self.es.part) // 2]
            for p in c.particles():
                rg += self.es.distance(p, com_particle)**2
            rg /= len(self.es.part)
            rg = np.sqrt(rg)
            self.assertAlmostEqual(c.radius_of_gyration(), rg, delta=1E-6)

            # Fractal dimension calc require gsl
            if not espressomd.has_features("GSL"):
                print("Skipping fractal dimension tests due to missing GSL dependency")
                return
            # The fractal dimension of a line should be 1

            dr = 0.
            self.assertAlmostEqual(
                c.fractal_dimension(dr=dr)[0], 1, delta=0.05)

            # Fractal dimension of a disk should be close to 2
            self.es.part.clear()
            center = np.array((0.1, .02, 0.15))
            for _ in range(3000):
                r_inv, phi = np.random.random(2) * np.array((0.2, 2 * np.pi))
                r = 1 / r_inv
                self.es.part.add(
                    pos=center + r * np.array((np.sin(phi), np.cos(phi), 0)))
            self.cs.clear()
            self.cs.run_for_all_pairs()
            cid = self.cs.cluster_ids()[0]
            df = self.cs.clusters[cid].fractal_dimension(dr=0.001)
            self.assertAlmostEqual(df[0], 2, delta=0.08)

    def test_analysis_for_bonded_particles(self):
        # Run cluster analysis
        self.cs.set_params(pair_criterion=BondCriterion(bond_type=0))
        self.cs.run_for_bonded_particles()

        # There should be one cluster containing particles 0 and 1
        self.assertEqual(len(self.cs.clusters), 1)
        self.assertEqual(
            self.cs.clusters[self.cs.cluster_ids()[0]].particle_ids(), [0, 1])

        # Check particle to cluster id mapping, once by ParticleHandle, once by
        # id
        self.assertEqual(self.cs.cid_for_particle(
            self.es.part[0]), self.cs.cluster_ids()[0])
        self.assertEqual(self.cs.cid_for_particle(1), self.cs.cluster_ids()[0])


if __name__ == "__main__":
    ut.main()
