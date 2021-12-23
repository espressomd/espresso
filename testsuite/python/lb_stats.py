# Copyright (C) 2010-2019 The ESPResSo project
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
              'viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5}

    system.periodicity = [1, 1, 1]
    system.time_step = 0.01
    system.cell_system.skin = 1.0
    dof = 3.

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def test_mass_momentum_thermostat(self):
        self.n_col_part = 100
        partcls = self.system.part.add(pos=np.random.random(
            (self.n_col_part, 3)) * self.system.box_l[0])
        if espressomd.has_features("MASS"):
            partcls.mass = 0.1 + np.random.random(
                len(self.system.part))

        self.system.thermostat.turn_off()

        self.lbf = self.lb_class(
            kT=self.params['temp'],
            visc=self.params['viscosity'],
            dens=self.params['dens'],
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
            for lb_node in self.lbf.nodes():
                dens = lb_node.density
                fluid_mass += dens
                fluid_temp += np.sum(lb_node.velocity**2) * dens

            # Normalize
            fluid_mass /= np.product(self.lbf.shape)
            fluid_temp *= self.system.volume() / (
                3. * np.product(self.lbf.shape)**2)

            # check mass conversation
            self.assertAlmostEqual(fluid_mass, self.params["dens"],
                                   delta=self.params["mass_prec_per_node"])

            # check momentum conservation
            # NOTE: this particle momentum prediction is due to the missing f/2 part in the
            #       LB fluid.
            particle_momentum = np.sum(
                [p.mass * p.v + 0.5 * p.f * self.system.time_step for p in self.system.part], axis=0)
            fluid_momentum = self.system.analysis.linear_momentum(False, True)
            np.testing.assert_allclose(
                particle_momentum + fluid_momentum, self.tot_mom,
                atol=self.params['mom_prec'])

            # Calc particle temperature
            e = self.system.analysis.energy()
            temp_particle = 2.0 / self.dof * e["kinetic"] / self.n_col_part

            # Update lists
            all_temp_particle.append(temp_particle)
            all_temp_fluid.append(fluid_temp)

        # import scipy.stats
        # temp_prec_particle = scipy.stats.norm.interval(0.95, loc=self.params["temp"],
        #   scale=np.std(all_temp_particle,ddof=1))[1] - self.params["temp"]
        # temp_prec_fluid = scipy.stats.norm.interval(0.95, loc=self.params["temp"],
        #   scale=np.std(all_temp_fluid,ddof=1))[1] -self.params["temp"]
        temp_prec_particle = 0.06 * self.params["temp"]
        temp_prec_fluid = 0.05 * self.params["temp"]

        self.assertAlmostEqual(
            np.mean(all_temp_fluid), self.params["temp"], delta=temp_prec_fluid)
        self.assertAlmostEqual(
            np.mean(all_temp_particle), self.params["temp"], delta=temp_prec_particle)


class TestLBCPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid
        self.params.update({"mom_prec": 1E-9, "mass_prec_per_node": 5E-8})


@utx.skipIfMissingGPU()
class TestLBGPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU
        self.params.update({"mom_prec": 1E-3, "mass_prec_per_node": 1E-5})


if __name__ == "__main__":
    ut.main()
