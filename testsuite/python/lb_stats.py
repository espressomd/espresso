#
# Copyright (C) 2010-2022 The ESPResSo project
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
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb


class TestLB:

    """
    Test the lattice-Boltzmann mass and momentum conservation.

    """
    system = espressomd.System(box_l=3 * [6.0])
    np.random.seed(1)
    params = {'tau': 0.01,
              'agrid': 0.5,
              'dens': 0.85,
              'kinematic_viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5}

    system.periodicity = [True, True, True]
    system.time_step = 0.01
    system.cell_system.skin = 0
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def test_mass_momentum_thermostat(self):
        self.n_col_part = 100
        # different types needed for lb_stats_hybrid test
        particles = self.system.part.add(
            type=self.n_col_part // 2 * [0, 1], pos=np.random.random(
                (self.n_col_part, 3)) * self.system.box_l[0])
        if espressomd.has_features("MASS"):
            particles.mass = 0.1 + np.random.random(
                len(self.system.part))

        self.system.thermostat.turn_off()

        self.lbf = self.lb_class(
            kT=self.params['temp'],
            kinematic_viscosity=self.params['kinematic_viscosity'],
            density=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0], seed=4)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            seed=3,
            gamma=self.params['friction'])
        # give particles a push
        for p in self.system.part:
            p.v = p.v + [0.1, 0.0, 0.0]

        self.fluidmass = self.params['dens']
        self.tot_mom = [0.0, 0.0, 0.0]
        for p in self.system.part:
            self.tot_mom += p.v * p.mass

        self.system.integrator.run(100)

        self.max_dmass = 0.0
        self.max_dm = [0, 0, 0]
        all_temp_particle = []
        all_temp_fluid = []

        # Integration
        for _ in range(20):
            self.system.integrator.run(15)

            # Summation vars
            fluid_mass = 0.0
            fluid_temp = 0.0

            # Go over lb lattice
            nodes_dens = self.lbf[:, :, :].density
            nodes_vel = np.sum(np.square(self.lbf[:, :, :].velocity), axis=3)
            fluid_mass += np.sum(nodes_dens)
            fluid_temp += np.sum(np.multiply(nodes_dens, nodes_vel))

            # Normalize
            fluid_mass /= np.product(self.lbf.shape)
            fluid_temp *= self.system.volume() / (
                3. * np.product(self.lbf.shape)**2)

            # check mass conversation
            self.assertAlmostEqual(fluid_mass, self.params["dens"], delta=1E-9)

            # check momentum conservation
            momentum = self.system.analysis.linear_momentum()
            f_2_correction = np.sum(
                self.system.part.all().f,
                axis=0) * self.system.time_step

            np.testing.assert_allclose(momentum + f_2_correction, self.tot_mom,
                                       atol=1E-10)

            temp_particle = np.average(
                [np.average(p.mass * p.v**2) for p in self.system.part])

            # Update lists
            all_temp_particle.append(temp_particle)
            all_temp_fluid.append(fluid_temp)

        # import scipy.stats
        # temp_prec_particle = scipy.stats.norm.interval(0.95, loc=self.params["temp"],
        #   scale=np.std(all_temp_particle,ddof=1))[1] - self.params["temp"]
        # temp_prec_fluid = scipy.stats.norm.interval(0.95, loc=self.params["temp"],
        #   scale=np.std(all_temp_fluid,ddof=1))[1] -self.params["temp"]
        temp_prec_particle = 0.05 * self.params["temp"]
        temp_prec_fluid = 0.05 * self.params["temp"]

        self.assertAlmostEqual(
            np.mean(all_temp_fluid), self.params["temp"], delta=temp_prec_fluid)
        self.assertAlmostEqual(
            np.mean(all_temp_particle), self.params["temp"], delta=temp_prec_particle)


@ut.skipIf(TestLB.n_nodes == 1,
           "LB with regular decomposition already tested with 2 MPI ranks")
@utx.skipIfMissingFeatures("WALBERLA")
class TestRegularLBWalberla(TestLB, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla

    def setUp(self):
        self.system.cell_system.set_regular_decomposition()


@utx.skipIfMissingFeatures("WALBERLA")
class TestNSquareLBWalberla(TestLB, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla

    def setUp(self):
        self.system.cell_system.set_n_square()


@utx.skipIfMissingFeatures("WALBERLA")
class TestHybrid0LBWalberla(TestLB, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla

    def setUp(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=0)
        self.params.update({"mom_prec": 1E-9, "mass_prec_per_node": 5E-8})


@utx.skipIfMissingFeatures("WALBERLA")
class TestHybrid1LBWalberla(TestLB, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla

    def setUp(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=0)
        self.params.update({"mom_prec": 1E-3, "mass_prec_per_node": 1E-5})


if __name__ == "__main__":
    ut.main()
