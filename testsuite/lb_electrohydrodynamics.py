import espressomd
import unittest as ut
import numpy as np

@ut.skipIf(not espressomd.has_features(["LB_ELECTROHYDRODYNAMICS"]),
           "Features not available, skipping test!")
class LBEHTest(ut.TestCase):
    from espressomd import lb
    s = espressomd.System(box_l=[6.0, 6.0, 6.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.params = {'time_step': 0.005,
                       'tau': 0.02,
                       'agrid': 0.5,
                       'dens': 0.85,
                       'viscosity': 30.0,
                       'friction': 3.0,
                       'temp': 0.0,
                       'skin': 0.2,
                       'muE': [0.1, 0.2, 0.3]}

        self.s.periodicity = [1, 1, 1]
        self.s.time_step = self.params['time_step']
        self.s.cell_system.skin = self.params['skin']

        for i in self.s.actors:
            self.s.actors.remove(i)

        self.lbf = self.lb.LBFluid(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.s.time_step,
            fric=self.params['friction'])

        self.s.actors.add(self.lbf)
        self.s.thermostat.set_lb(kT=self.params['temp'])

    def test(self):
        s = self.s

        s.part.add(pos=0.5 * s.box_l, mu_E=self.params['muE'])

        mu_E = np.array(self.params['muE'])
        # Terminal velocity is mu_E minus the momentum the fluid
        # got by accelerating the particle in the beginning.
        v_term =  (1. - 1. / (s.box_l[0]*s.box_l[1]*s.box_l[2]*self.params['dens']))*mu_E

        s.integrator.run(steps=1000)

        np.testing.assert_allclose(v_term, np.copy(s.part[0].v), atol=1e-5)

if __name__ == "__main__":
    ut.main()

