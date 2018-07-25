from __future__ import print_function

import itertools
import unittest as ut
import numpy as np

import espressomd
import espressomd.lb
from espressomd import *
from tests_common import abspath


class TestLB(ut.TestCase):

    """
    Basic tests of the Lattice Boltzmann implementation

    * mass and momentum conservation
    * temperature
    * particle viscous coupling
    * application of external force densities
    * setting and retrieving lb node velocities

    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = range(n_nodes)
    params = {'int_steps': 25,
              'int_times': 10,
              'time_step': 0.01,
              'tau': 0.02,
              'agrid': 0.5,
              'box_l': 12.0,
              'dens': 0.85,
              'viscosity': 30.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5,
              'skin': 0.2,
              'temp_confidence': 10}
    if espressomd.has_features("SHANCHEN"):
       params.update({"dens":2*[params["dens"]]})
     
    if espressomd.has_features("ROTATION"):
        dof = 6.
    else:
        dof = 3.

    system.box_l = [
        params['box_l'],
        params['box_l'],
        params['box_l']]
    system.periodicity = [1, 1, 1]
    system.time_step = params['time_step']
    system.cell_system.skin = params['skin']

    def test_mass_momentum_thermostat(self):
        self.system.actors.clear()
        self.system.part.clear()
        # import particle data
        self.data = np.genfromtxt(abspath("data/lb_system.data"))

        for particle in self.data:
            id = particle[0]
            typ = particle[1]
            pos = particle[3:6]
            f = particle[9:]
            v = particle[6:9]
            p = self.system.part.add(id=int(id), pos=pos, v=v, type=int(typ))
            if espressomd.has_features("ROTATION"):
                p.rotation = [1, 1, 1]

        self.n_col_part = len(self.system.part)

        self.system.thermostat.set_langevin(
            kT=self.params['temp'], gamma=self.params['gamma'])
        self.system.integrator.run(50)
        # kill particle motion
        for i in range(self.n_col_part):
            self.system.part[i].v = [0.0, 0.0, 0.0]
        self.system.thermostat.turn_off()

        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            fric=self.params['friction'], ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(kT=self.params['temp'])
        # give particles a push
        for i in range(self.n_col_part):
            self.system.part[i].v = self.system.part[i].v + [0.1, 0.0, 0.0]

        self.fluidmass = self.params['dens']
        self.tot_mom = [0.0, 0.0, 0.0]
        for i in range(self.n_col_part):
            self.tot_mom = self.tot_mom + self.system.part[i].v

        self.system.integrator.run(50)

        self.max_dmass = 0.0
        self.max_dm = [0, 0, 0]
        self.avg_temp = 0.0
        self.avg_fluid_temp = 0.0

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

            # check momentum conservation
            c_mom = self.system.analysis.analyze_linear_momentum()
            dm = abs(c_mom - self.tot_mom)
            for j in range(3):
                if dm[j] > self.max_dm[j]:
                    self.max_dm[j] = dm[j]
            self.assertTrue(
                self.max_dm[0] <= self.params['mom_prec'] and self.max_dm[1] <= self.params[
                    'mom_prec'] and self.max_dm[2] <= self.params['mom_prec'],
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

    def test_set_get_u(self):
        self.system.actors.clear()
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            fric=self.params['friction'], ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        v_fluid = np.array([1.2, 4.3, 0.2])
        self.lbf[0, 0, 0].velocity = v_fluid
        np.testing.assert_allclose(
            np.copy(self.lbf[0, 0, 0].velocity), v_fluid, atol=1e-4)

    def test_viscous_coupling(self):
        self.system.thermostat.turn_off()
        self.system.actors.clear()
        self.system.part.clear()
        v_part = np.array([1, 2, 3])
        v_fluid = np.array([1.2, 4.3, 0.2])
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            fric=self.params['friction'], ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        self.system.part.add(
            pos=[0.5 * self.params['agrid']] * 3, v=v_part, fix=[1, 1, 1])
        self.lbf[0, 0, 0].velocity = v_fluid
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(self.system.part[0].f), -self.params['friction'] * (v_part - v_fluid),atol=1E-6)

    def test_a_ext_force_density(self):
        self.system.thermostat.turn_off()
        self.system.actors.clear()
        self.system.part.clear()
        ext_force_density = [2.3, 1.2, 0.1]
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            fric=self.params['friction'], ext_force_density=ext_force_density)
        self.system.actors.add(self.lbf)
        n_time_steps = 5
        self.system.integrator.run(n_time_steps)
        # ext_force_density is a force density, therefore v = ext_force_density / dens * tau * (n_time_steps - 0.5)
        # (force is applied only to the second half of the first integration step)
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps - 0.5) / self.params['dens']
        for n in list(itertools.combinations(range(int(self.system.box_l[0] / self.params['agrid'])), 3)):
            np.testing.assert_allclose(
                np.copy(self.lbf[n].velocity), fluid_velocity,atol=1E-6)

