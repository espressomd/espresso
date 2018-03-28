# Basic tests of the Lattice Boltzmann implementation
#
# 1) check conservation of fluid mass
# 2) check conservation of total momentum
# 3) measure temperature of colloid and fluid

from __future__ import print_function
import sys
import numpy as np
import unittest as ut
import math

import espressomd
import espressomd.lb
import espressomd.shapes
from tests_common import abspath


@ut.skipIf(not espressomd.has_features(["LB_GPU", "LENNARD_JONES"]) or espressomd.has_features(
    "SHANCHEN"), "Features not available, skipping test!")
class TestLBGPU(ut.TestCase):
    box_l = 10.0
    system = espressomd.System(box_l=[box_l] * 3)
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = range(n_nodes)
    int_steps = 5
    int_times = 10
    time_step = 0.005
    tau = 0.02
    agrid = 1.0
    dens = 0.85
    viscosity = 30.0
    friction = 2.0
    temp = 1.0
    gamma = 1.0
    skin = 0.4
    mom_prec = 1.e-2
    mass_prec_per_node = 4.e-8
    temp_confidence = 10

    def setUp(self):
        # setup parameters
        system = self.system
        system.part.clear()
        system.time_step = self.time_step

        if espressomd.has_features("ROTATION"):
            self.dof = 6.
        else:
            self.dof = 3.

        system.periodicity = [1, 1, 1]
        system.time_step = self.time_step
        system.cell_system.skin = self.skin

    def tearDown(self):
        self.system.actors.clear()

    def setup_lb_fluid(self, ext_force =[0.0, 0.0, 0.0]):
        lbf = espressomd.lb.LBFluidGPU(
            visc=self.viscosity,
            dens=self.dens,
            agrid=self.agrid,
            tau=self.system.time_step,
            fric=self.friction,
            ext_force = ext_force)
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

    def test_thermostat(self):
        system = self.system
        self.place_particles()
        system.thermostat.set_langevin(kT=self.temp, gamma=self.gamma)
        system.integrator.run(100)
        # kill particle motion
        n_col_part = len(system.part)
        for i in range(n_col_part):
            system.part[i].v = [0.0, 0.0, 0.0]
        system.thermostat.turn_off()
        lbf = self.setup_lb_fluid()
        system.thermostat.set_lb(kT=self.temp)
        # give particles a push
        for i in range(n_col_part):
            system.part[i].v = system.part[i].v + [0.1, 0.0, 0.0]

        fluidmass = self.dens
        tot_mom = [0.0, 0.0, 0.0]
        for i in range(n_col_part):
            tot_mom += system.part[i].v

        system.integrator.run(100)

        max_dmass = 0.0
        max_dm = [0, 0, 0]
        avg_temp = 0.0
        avg_fluid_temp = 0.0

        # Integration
        for i in range(self.int_times):
            system.integrator.run(self.int_steps)
            fluidmass_sim = 0.0
            fluid_temp_sim = 0.0
            node_v_list = []
            node_dens_list = []
            for i in range(int(self.box_l / self.agrid)):
                for j in range(int(self.box_l / self.agrid)):
                    for k in range(int(self.box_l / self.agrid)):
                        node_v_list.append(lbf[i, j, k].velocity)
                        node_dens_list.append(lbf[i, j, k].density[0])

            # check mass conservation
            fluidmass_sim = sum(node_dens_list) / len(node_dens_list)

            dmass = abs(fluidmass_sim - fluidmass)
            if dmass > max_dmass:
                max_dmass = dmass
            self.assertTrue(
                max_dmass < self.mass_prec_per_node,
                msg="Fluid mass deviation too high\ndeviation: {} accepted deviation: {}".format(
                    max_dmass,
                    self.mass_prec_per_node))

            # check momentum conservation
            c_mom = system.analysis.analyze_linear_momentum()
            dm = abs(c_mom - tot_mom)
            for j in range(3):
                if dm[j] > max_dm[j]:
                    max_dm[j] = dm[j]
            self.assertTrue(
                max_dm[0] <= self.mom_prec and max_dm[1] <= self.mom_prec and max_dm[2] <= self.mom_prec,
                msg="Momentum deviation too high\ndeviation: {} accepted deviation: {}".format(
                    max_dm,
                    self.mom_prec))

            # check temp of particles
            e = system.analysis.energy()
            temp_particle = 2.0 / self.dof * e["total"] / n_col_part
            avg_temp += temp_particle / self.int_times
            # check temp of fluid
            fluid_temp = 0
            for j in range(len(node_dens_list)):
                fluid_temp += (1.0 / 3.0) * (node_v_list[j][0]**2.0 + node_v_list[j][1]**2.0 +
                                             node_v_list[j][2]**2.0) * node_dens_list[j] * (self.box_l)**3 / len(node_dens_list)**2
            avg_fluid_temp += fluid_temp / self.int_times

        temp_dev = (2.0 / (n_col_part * 3.0))**0.5
        temp_prec = self.temp_confidence * temp_dev / (self.int_times)**0.5

        self.assertTrue(
            abs(
                avg_temp -
                self.temp) < temp_prec,
            msg="Particle temperature deviation too high\ndeviation: {} accepted deviation: {}".format(
                abs(
                    avg_temp -
                    self.temp),
                temp_prec))
        self.assertTrue(
            abs(
                avg_fluid_temp -
                self.temp) < temp_prec,
            msg="Fluid temperature deviation too high\ndeviation: {} accepted deviation: {}".format(
                abs(
                    avg_fluid_temp -
                    self.temp),
                temp_prec))

    def test_zboundary_slip(self):
        lbf = self.setup_lb_fluid(ext_force = [0.0, 0.0, 0.1])
        wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=1.0)
        wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-9.0)
        boundary1 = espressomd.lbboundaries.LBBoundary(shape=wall_shape1)
        boundary2 = espressomd.lbboundaries.LBBoundary(shape=wall_shape2)
        self.system.lbboundaries.add(boundary1)
        self.system.lbboundaries.add(boundary2)
        self.system.integrator.run(100)
        zero_slip_positions = []
        for y in range(int(math.floor(self.system.box_l[1]/self.agrid))):
            zero_slip_positions.append([self.agrid, (y+0.5)*self.agrid, self.agrid])
        velocities = lbf.get_interpolated_fluid_velocity_at_positions(np.array(zero_slip_positions))
        np.testing.assert_array_equal(velocities, np.zeros_like(velocities))



if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(TestLBGPU))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
