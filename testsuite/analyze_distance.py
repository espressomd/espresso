from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd


@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeDistance(ut.TestCase):
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(1234)

    @classmethod
    def setUpClass(self):
        box_l = 50.0
        self.system.box_l = [box_l, box_l, box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0,
            cutoff=2**(1. / 6.), shift="auto")
        self.system.thermostat.set_langevin(kT=1., gamma=1.)
        for i in range(100):
            self.system.part.add(id=i, pos=np.random.random(3) * box_l)
        self.system.force_cap = 10
        i = 0
        min_dist = self.system.analysis.min_dist()
        while (i < 50 and min_dist < 0.9):
            system.integrator.run(100)
            min_dist = self.system.analysis.min_dist()
            i += 1
            lj_cap = lj_cap + 10
            self.system.force_cap = lj_cap
        self.system.force_cap = 0
        self.system.integrator.run(1000)

    # python version of the espresso core function
    def min_dist(self):
        r = np.array(self.system.part[:].pos)
        # this generates indices for all i<j combinations
        ij = np.triu_indices(len(r), k=1)
        r_ij = np.fabs(r[ij[0]] - r[ij[1]])
        # check smaller distances via PBC
        r_ij = np.where(r_ij > 0.5 * self.system.box_l, self.system.box_l-r_ij, r_ij)
        dist = np.sum(r_ij**2, axis=1)
        return np.sqrt(np.min(dist))

    # python version of the espresso core function
    def nbhood(self, pos, r_catch):
        dist = np.fabs(np.array(self.system.part[:].pos) - pos)
        # check smaller distances via PBC
        dist = np.where(dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=1)
        return np.where(dist < r_catch**2)[0]

    # python version of the espresso core function, using pos
    def dist_to_pos(self, pos):
        dist = np.fabs(self.system.part[:].pos - pos)
        # check smaller distances via PBC
        dist = np.where(dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=-1)
        return np.sqrt(np.min(dist))

    # python version of the espresso core function, using id
    def dist_to_id(self, id):
        dist = np.fabs(np.delete(self.system.part[:].pos, id, axis=0) - self.system.part[id].pos)
        # check smaller distances via PBC
        dist = np.where(dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=1)
        return np.sqrt(np.min(dist))

    def test_min_dist(self):
        # try five times
        for i in range(5):
            self.assertAlmostEqual(self.system.analysis.min_dist(),
                                   self.min_dist(),
                                   delta=1e-7)
            self.system.integrator.run(100)

    def test_nbhood(self):
        # try five times
        for i in range(1, 10, 2):
            self.assertTrue(
                np.allclose(self.system.analysis.nbhood([i, i, i], i * 2),
                            self.nbhood([i, i, i], i * 2)))
            self.system.integrator.run(100)

    def test_dist_to_pos(self):
        # try five times
        for i in range(5):
            self.assertTrue(
                np.allclose(self.system.analysis.dist_to(pos=[i, i, i]),
                            self.dist_to_pos([i, i, i])))
            self.system.integrator.run(100)

    def test_dist_to_id(self):
        # try five times
        for i in range(5):
            self.assertAlmostEqual(self.system.analysis.dist_to(id=i),
                                   self.dist_to_id(i))
            self.system.integrator.run(100)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
