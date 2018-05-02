from __future__ import print_function

import sys
import unittest as ut
import numpy as np
import math

import espressomd
import espressomd.lb
import espressomd.shapes
from tests_common import abspath


@ut.skipIf(not espressomd.has_features(["LB"]),
           "Features not available, skipping test!")
class LBTest(ut.TestCase):
    """
    Basic tests of the Lattice Boltzmann implementation

    1) check conservation of fluid mass
    2) check conservation of total momentum
    3) measure temperature of colloid and fluid
    4) no-slip boundary condition

    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = range(n_nodes)
    params = {'int_steps': 50,
              'int_times': 10,
              'time_step': 0.01,
              'tau': 0.02,
              'agrid': 1.0,
              'box_l': 10.0,
              'dens': 0.85,
              'viscosity': 30.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5,
              'skin': 0.2,
              'mom_prec_gpu': 1.e-3,
              'mom_prec_cpu': 1.e-11,
              'mass_prec_per_node': 4.e-8,
              'temp_confidence': 10}

    def setUp(self):
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

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def setup_lb_fluid(self, gpu_LB=False, ext_force = [0.0, 0.0, 0.0]):
        if gpu_LB:
            lbf = espressomd.lb.LBFluidGPU(
                visc=self.params['viscosity'],
                dens=self.params['dens'],
                agrid=self.params['agrid'],
                tau=self.system.time_step,
                fric=self.params['friction'],
                ext_force=ext_force)
            self.system.actors.add(lbf)
        else:
            lbf = espressomd.lb.LBFluid(
                visc=self.params['viscosity'],
                dens=self.params['dens'],
                agrid=self.params['agrid'],
                tau=self.system.time_step,
                fric=self.params['friction'],
                ext_force=ext_force)
            self.system.actors.add(lbf)
        return lbf

    def place_particles(self):
        data = np.genfromtxt(abspath("data/lb_system.data"))
        for particle in data:
            id = particle[0]
            typ = particle[1]
            pos = particle[3:6]
            f = particle[9:]
            v = particle[6:9]
            p = self.system.part.add(id=int(id), pos=pos, v=v, type=int(typ))
            if espressomd.has_features("ROTATION"):
                p.rotation = [1, 1, 1]

    def perform_thermostat_test(self, gpu_LB=False):
        system = self.system
        self.place_particles()
        system.thermostat.set_langevin(kT=self.params['temp'], gamma=self.params['gamma'])
        system.integrator.run(100)
        # kill particle motion
        n_col_part = len(system.part)
        for i in range(n_col_part):
            system.part[i].v = [0.0, 0.0, 0.0]
        system.thermostat.turn_off()
        lbf = self.setup_lb_fluid(gpu_LB=gpu_LB)
        system.thermostat.set_lb(kT=self.params['temp'])
        # give particles a push
        for i in range(n_col_part):
            system.part[i].v = system.part[i].v + [0.1, 0.0, 0.0]

        fluidmass = self.params['dens']
        tot_mom = [0.0, 0.0, 0.0]
        for i in range(n_col_part):
            tot_mom += system.part[i].v

        system.integrator.run(100)

        max_dmass = 0.0
        max_dm = [0, 0, 0]
        avg_temp = 0.0
        avg_fluid_temp = 0.0

        # Integration
        for i in range(self.params['int_times']):
            system.integrator.run(self.params['int_steps'])
            fluidmass_sim = 0.0
            fluid_temp_sim = 0.0
            node_v_list = []
            node_dens_list = []
            for i in range(int(self.params['box_l'] / self.params['agrid'])):
                for j in range(int(self.params['box_l'] / self.params['agrid'])):
                    for k in range(int(self.params['box_l'] / self.params['agrid'])):
                        node_v_list.append(lbf[i, j, k].velocity)
                        node_dens_list.append(lbf[i, j, k].density[0])

            # check mass conservation
            fluidmass_sim = sum(node_dens_list) / len(node_dens_list)

            dmass = abs(fluidmass_sim - fluidmass)
            if dmass > max_dmass:
                max_dmass = dmass
            self.assertTrue(
                max_dmass < self.params['mass_prec_per_node'],
                msg="Fluid mass deviation too high\ndeviation: {} accepted deviation: {}".format(
                    max_dmass,
                    self.params['mass_prec_per_node']))

            # check momentum conservation
            c_mom = system.analysis.analyze_linear_momentum()
            dm = abs(c_mom - tot_mom)
            for j in range(3):
                if dm[j] > max_dm[j]:
                    max_dm[j] = dm[j]
            if gpu_LB:
                mom_prec = self.params['mom_prec_gpu']
            else:
                mom_prec = self.params['mom_prec_cpu']
            self.assertTrue(
                max_dm[0] <= mom_prec and max_dm[1] <= mom_prec and max_dm[2] <= mom_prec,
                msg="Momentum deviation too high\ndeviation: {} accepted deviation: {}".format(
                    max_dm,
                    mom_prec))

            # check temp of particles
            e = system.analysis.energy()
            temp_particle = 2.0 / self.dof * e["total"] / n_col_part
            avg_temp += temp_particle / self.params['int_times']
            # check temp of fluid
            fluid_temp = 0
            for j in range(len(node_dens_list)):
                fluid_temp += (1.0 / 3.0) * (node_v_list[j][0]**2.0 + node_v_list[j][1]**2.0 +
                                             node_v_list[j][2]**2.0) * node_dens_list[j] * (self.params['box_l'])**3 / len(node_dens_list)**2
            avg_fluid_temp += fluid_temp / self.params['int_times']

        temp_dev = (2.0 / (n_col_part * 3.0))**0.5
        temp_prec = self.params['temp_confidence'] * temp_dev / (self.params['int_times'])**0.5

        self.assertTrue(
            abs(
                avg_temp -
                self.params['temp']) < temp_prec,
            msg="Particle temperature deviation too high\ndeviation: {} accepted deviation: {}".format(
                abs(
                    avg_temp -
                    self.params['temp']),
                temp_prec))
        self.assertTrue(
            abs(
                avg_fluid_temp -
                self.params['temp']) < temp_prec,
            msg="Fluid temperature deviation too high\ndeviation: {} accepted deviation: {}".format(
                abs(
                    avg_fluid_temp -
                    self.params['temp']),
                temp_prec))

    def test_thermostat_cpu(self):
        self.perform_thermostat_test()
    
    @ut.skipIf(not espressomd.has_features("LB_GPU"), "Skipping test because LB_GPU feature not compiled in.")
    def test_thermostat_gpu(self):
        self.perform_thermostat_test(gpu_LB=True)


    def perform_zboundary_slip_test(self, gpu_LB=False):
        lbf = self.setup_lb_fluid(ext_force = [0.0, 0.0, 0.1], gpu_LB=gpu_LB)
        wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=1.0)
        wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-9.0)
        boundary1 = espressomd.lbboundaries.LBBoundary(shape=wall_shape1)
        boundary2 = espressomd.lbboundaries.LBBoundary(shape=wall_shape2)
        self.system.lbboundaries.add(boundary1)
        self.system.lbboundaries.add(boundary2)
        self.system.integrator.run(100)
        zero_slip_positions = []
        for y in range(int(math.floor(self.system.box_l[1] / self.params['agrid']))):
            zero_slip_positions.append([self.params['agrid'], (y+0.5)*self.params['agrid'], self.params['agrid']])
        if gpu_LB:
            velocities = lbf.get_interpolated_fluid_velocity_at_positions(np.array(zero_slip_positions))
        else:
            velocities = [lbf.get_interpolated_velocity_local(p) for p in zero_slip_positions]
        np.testing.assert_array_equal(velocities, np.zeros_like(velocities))

    def test_zboundary_slip_cpu(self):
        self.perform_zboundary_slip_test()

    @ut.skipIf(not espressomd.has_features("LB_GPU"), "Skipping test because LB_GPU feature not compiled in.")
    def test_zboundary_slip_gpu(self):
        self.perform_zboundary_slip_test(gpu_LB=True)

if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(LBTest))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
