#
# Copyright (C) 2019 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import espressomd.observables


class Observables(ut.TestCase):
    n_tries = 50
    n_parts = 8
    box_l = 5.
    system = espressomd.System(box_l=3 * [box_l])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

    @classmethod
    def setUpClass(cls):
        for i in range(cls.n_parts):
            cls.system.part.add(pos=[1 + i, 1 + i, 1 + i], id=i)

    def scatter_particles(self):
        pos = np.random.uniform(low=-1.5 * self.box_l, high=1.5 * self.box_l,
                                size=[self.n_parts, 3])
        for i in range(self.n_parts):
            self.system.part[i].pos = pos[i]
        return pos

    def test_ParticleDistances(self):
        """
        Check ParticleDistances for a single particle pair and a chain.
        """
        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.ParticleDistances(ids=[0, 1])
        obs_chain = espressomd.observables.ParticleDistances(ids=pids)
        for _ in range(self.n_tries):
            pos = self.scatter_particles()
            distances = []
            for i in pids[:-1]:
                distances.append(np.linalg.norm(pos[i + 1] - pos[i]))
            res_obs_single = obs_single.calculate()
            res_obs_chain = obs_chain.calculate()
            self.assertEqual(obs_single.n_values(), 1)
            self.assertEqual(obs_chain.n_values(), self.n_parts - 1)
            self.assertAlmostEqual(res_obs_single[0], distances[0], places=9)
            np.testing.assert_array_almost_equal(
                res_obs_chain, distances, decimal=9,
                err_msg="Data did not agree for observable ParticleDistances")

    def test_ParticleAngles(self):
        """
        Check ParticleAngles for a single particle triple and a chain.
        """
        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.ParticleAngles(ids=[0, 1, 2])
        obs_chain = espressomd.observables.ParticleAngles(ids=pids)
        for _ in range(self.n_tries):
            pos = self.scatter_particles()
            vectors = pos[1:] - pos[:-1]
            angles = []
            for i in pids[:-2]:
                v1 = vectors[i]
                v2 = vectors[i + 1]
                angles.append(np.arccos(
                    -np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
            res_obs_single = obs_single.calculate()
            res_obs_chain = obs_chain.calculate()
            self.assertEqual(obs_single.n_values(), 1)
            self.assertEqual(obs_chain.n_values(), self.n_parts - 2)
            self.assertAlmostEqual(res_obs_single[0], angles[0], places=9)
            np.testing.assert_array_almost_equal(
                res_obs_chain, angles, decimal=9,
                err_msg="Data did not agree for observable ParticleAngles")

    def test_ParticleDihedrals(self):
        """
        Check ParticleDihedrals, for a single particle quadruple and a chain.
        """
        pids = list(range(self.n_parts))
        obs_single = espressomd.observables.ParticleDihedrals(ids=[0, 1, 2, 3])
        obs_chain = espressomd.observables.ParticleDihedrals(ids=pids)
        counter = 0
        loops = 0
        while counter < self.n_tries:
            loops += 1
            self.assertTrue(loops < self.n_tries + 5, "np.random.uniform() "
                            "placed particles such that dihedral angles were undefined")
            pos = self.scatter_particles()
            vectors = pos[1:] - pos[:-1]
            dihedrals = []
            for i in pids[:-3]:
                v1 = vectors[i]
                v2 = vectors[i + 1]
                v3 = vectors[i + 2]
                c1 = np.cross(v1, v2)
                c2 = np.cross(v2, v3)
                # only carry out test if dihedral angle is defined
                if np.linalg.norm(c1) < 1e-8 or np.linalg.norm(c2) < 1e-8:
                    continue
                dihedrals.append(np.arctan2(
                    np.dot(np.cross(c1, c2), v2) / np.linalg.norm(v2),
                    np.dot(c1, c2)))
            res_obs_single = obs_single.calculate()
            res_obs_chain = obs_chain.calculate()
            self.assertEqual(obs_single.n_values(), 1)
            self.assertEqual(obs_chain.n_values(), self.n_parts - 3)
            self.assertAlmostEqual(res_obs_single[0], dihedrals[0], places=9)
            np.testing.assert_array_almost_equal(
                res_obs_chain, dihedrals, decimal=9,
                err_msg="Data did not agree for observable ParticleDihedrals")
            counter += 1

if __name__ == "__main__":
    ut.main()
