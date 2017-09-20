
from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
import espressomd.lb
from espressomd import *
from tests_common import abspath


@ut.skipIf(not espressomd.has_features(["LB", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class LBTest(ut.TestCase):
    """
    Basic tests of the Lattice Boltzmann implementation
    
    1) check conservation of fluid mass
    2) check conservation of total momentum
    3) measure temperature of colloid and fluid

    """
    system = espressomd.System()
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = np.random.randint(low=1, high=2**31 - 1, size=n_nodes)

    def setUp(self):
        self.params = {'int_steps': 100,
                       'int_times': 20,
                       'time_step': 0.005,
                       'tau': 0.02,
                       'agrid': 1.0,
                       'box_l': 30.0,
                       'dens': 0.85,
                       'viscosity': 30.0,
                       'friction': 2.0,
                       'temp': 1.0,
                       'gamma': 1.0,
                       'skin': 0.4,
                       'mom_prec': 1.e-11,
                       'mass_prec_per_node': 4.e-8,
                       'temp_confidence': 10}

        if espressomd.has_features("ROTATION"):
            self.dof = 6.
        else:
            self.dof = 3.

        self.system.box_l = [
            self.params['box_l'],
            self.params['box_l'],
            self.params['box_l']]
        self.system.periodicity = [1, 1, 1]
        self.system.time_step = self.params['time_step']
        self.system.cell_system.skin = self.params['skin']

        # clear actors that might be left from prev tests
        for i in self.system.actors:
            self.system.actors.remove(i)
        self.system.part.clear()
        # import particle data
        self.data = np.genfromtxt(abspath("data/lb_system.data"))

        for particle in self.data:
            id = particle[0]
            typ = particle[1]
            pos = particle[3:6]
            f = particle[9:]
            v = particle[6:9]
            self.system.part.add(id=int(id), pos=pos, v=v, type=int(typ))

        self.n_col_part = len(self.system.part)

        self.system.thermostat.set_langevin(
            kT=self.params['temp'], gamma=self.params['gamma'])
        self.system.integrator.run(1000)
        # kill particle motion
        for i in range(self.n_col_part):
            self.system.part[i].v = [0.0, 0.0, 0.0]
        self.system.thermostat.turn_off()

        self.lbf = lb.LBFluid(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            fric=self.params['friction'])
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(kT=self.params['temp'])
        # give particles a push
        for i in range(self.n_col_part):
            self.system.part[i].v = self.system.part[i].v + [0.1, 0.0, 0.0]

        self.fluidmass = self.params['dens']
        self.tot_mom = [0.0, 0.0, 0.0]
        for i in range(self.n_col_part):
            self.tot_mom = self.tot_mom + self.system.part[i].v

        self.system.integrator.run(1000)

        self.max_dmass = 0.0
        self.max_dm = [0, 0, 0]
        self.avg_temp = 0.0
        self.avg_fluid_temp = 0.0

    def test(self):
        # Integration
        for i in range(self.params['int_times']):
            self.system.integrator.run(self.params['int_steps'])
            fluidmass_sim = 0.0
            fluid_temp_sim = 0.0
            node_v_list = []
            node_dens_list = []
            n_nodes = int(self.params['box_l'] / self.params['agrid'])
            for i in range(n_nodes):
                for j in range(n_nodes):
                    for k in range(n_nodes):
                        node_v_list.append(self.lbf[i, j, k].velocity)
                        node_dens_list.append(self.lbf[i, j, k].density[0])

            # check mass conversation
            fluidmass_sim = sum(node_dens_list) / len(node_dens_list)

            dmass = abs(fluidmass_sim - self.fluidmass)  # /len(node_dens_list)
            if dmass > self.max_dmass:
                self.max_dmass = dmass
            self.assertTrue(
                self.max_dmass < self.params['mass_prec_per_node'],
                msg="fluid mass deviation too high\ndeviation: {}   accepted deviation: {}".format(
                    self.max_dmass,
                    self.params['mass_prec_per_node']))

            # check momentum conversation
            c_mom = self.system.analysis.analyze_linear_momentum()
            dm = abs(c_mom - self.tot_mom)
            #self.assertTrue(dm[0] <= mom_prec and dm[1] <= mom_prec and dm[2] <= mom_prec)
            for j in range(3):
                if dm[j] > self.max_dm[j]:
                    self.max_dm[j] = dm[j]
            self.assertTrue(
                self.max_dm[0] <= self.params['mom_prec'] and self.max_dm[1] <= self.params['mom_prec'] and self.max_dm[2] <= self.params['mom_prec'],
                msg="momentum deviation too high\ndeviation: {}  accepted deviation: {}".format(
                    self.max_dm,
                    self.params['mom_prec']))

            # check temp of particles
            e = self.system.analysis.energy()
            temp_particle = 2.0 / self.dof * e["kinetic"] / self.n_col_part
            self.avg_temp = self.avg_temp + \
                temp_particle / self.params['int_times']
            # check temp of fluid
            fluid_temp = 0
            for j in range(len(node_dens_list)):
                fluid_temp = fluid_temp + (1.0 / 3.0) * (node_v_list[j][0]**2.0 + node_v_list[j][1] ** 2.0 +
                                                         node_v_list[j][2]**2.0) * node_dens_list[j] * (self.params['box_l'])**3 / len(node_dens_list)**2
            self.avg_fluid_temp = self.avg_fluid_temp + \
                fluid_temp / self.params['int_times']

        temp_dev = (2.0 / (self.n_col_part * 3.0))**0.5
        temp_prec = self.params['temp_confidence'] * \
            temp_dev / (self.params['int_times'])**0.5

        print("maximal mass deviation: {}  accepted deviation: {}".format(
            self.max_dmass, self.params['mass_prec_per_node']))
        print("maximal momentum deviation: {}  accepted deviation: {}".format(
            self.max_dm, self.params['mom_prec']))
        print("average temperature: {}".format(self.avg_temp))
        print("average fluid temperature: {}".format(self.avg_fluid_temp))
        print("set temperature: {}".format(self.params['temp']))
        print("maximally accepted deviation: {}".format(temp_prec))

        self.assertTrue(
            abs(
                self.avg_temp -
                self.params['temp']) < temp_prec,
            msg="particle temperature deviation too high\ndeviation: {}  accepted deviation: {}".format(
                abs(
                    self.avg_temp -
                    self.params['temp']),
                temp_prec))
        self.assertTrue(
            abs(
                self.avg_fluid_temp -
                self.params['temp']) < temp_prec,
            msg="fluid temperature deviation too high\ndeviation: {}  accepted deviation: {}".format(
                abs(
                    self.avg_fluid_temp -
                    self.params['temp']),
                temp_prec))


if __name__ == "__main__":
    ut.main()
