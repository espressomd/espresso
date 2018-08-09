from __future__ import print_function
import unittest as ut
import espressomd

@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]),
           "Features not available, skipping test!")
class PairTest(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.s.time_step = 0.1
        self.s.thermostat.turn_off()

        self.s.part.clear()
        self.s.box_l = 3*[10.]
        self.s.cell_system.skin = 0.3

        # Force an appropriate cell grid
        self.s.non_bonded_inter[0,0].lennard_jones.set_params(
                epsilon=0.0, sigma=1.0,
                cutoff=1.5, shift=0.0)

        vel = [1., 2., 3.]

        self.s.part.add(id=0, pos=[5., 4.5, 5.], v=vel)
        self.s.part.add(id=1, pos=[5., 5.5, 5.], v=vel)
        self.s.part.add(id=2, pos=[9.5, 5., 5.], v=vel)
        self.s.part.add(id=3, pos=[0.5, 5., 5.], v=vel)
        self.s.part.add(id=4, pos=[5., 5., 9.5], v=vel)
        self.s.part.add(id=5, pos=[5., 5., 0.5], v=vel)
        self.s.part.add(id=6, pos=[5., 9.5, 5.], v=vel)
        self.s.part.add(id=7, pos=[5., 0.5, 5.], v=vel)
        self.s.part.add(id=8, pos=[5., 9.5, 9.5], v=vel)
        self.s.part.add(id=9, pos=[5., 0.5, 0.5], v=vel)
        self.s.part.add(id=10, pos=[1., 1., 1.], v=vel)
        self.s.part.add(id=11, pos=[9., 9., 9.], v=vel)

    def expected_pairs(self,periodicity):
        if all(periodicity == (1,1,1)):
            return [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
        elif all(periodicity == (1,1,0)):
            return [(0, 1), (2, 3), (6, 7)]

    def check(self):
        pairs = self.s.cell_system.get_pairs_(1.5)
        epairs = self.expected_pairs(self.s.periodicity)
        print(pairs)

        self.assertEqual(len(pairs), len(epairs))
        for p in pairs:
            self.assertTrue(p in epairs)

    def test_nsquare(self):
        self.s.cell_system.set_n_square()
        self.s.periodicity = 1,1,1

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

    @ut.skipIf(not espressomd.has_features(["PARTIAL_PERIODIC"]),
       "Features not available, skipping test!")
    def test_nsquare_partial_z(self):
        self.s.cell_system.set_n_square()
        self.s.periodicity = 1,1,0

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

    def test_dd(self):
        self.s.cell_system.set_domain_decomposition()
        self.s.cell_system.min_num_cells = 8
        self.s.periodicity = 1,1,1

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

    @ut.skipIf(not espressomd.has_features(["PARTIAL_PERIODIC"]),
       "Features not available, skipping test!")
    def test_dd_partial_z(self):
        self.s.cell_system.set_domain_decomposition()
        self.s.periodicity = 1,1,0

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

    def test_layered(self):
        self.s.cell_system.set_layered()
        self.s.periodicity = 1,1,1

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

    @ut.skipIf(not espressomd.has_features(["PARTIAL_PERIODIC"]),
       "Features not available, skipping test!")
    def test_layered_partial_z(self):
        self.s.cell_system.set_layered()
        self.s.periodicity = 1,1,0

        self.s.integrator.run(0)
        self.check()
        self.s.integrator.run(100)
        self.check()

if __name__ == "__main__":
    ut.main()
