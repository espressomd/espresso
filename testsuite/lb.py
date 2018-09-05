# Copyright (C) 2010-2018 The ESPResSo project
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

import itertools
import unittest as ut
import numpy as np

import espressomd
import espressomd.lb
from tests_common import abspath


class TestLB(object):

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
    #np.random.seed = 1
    params = {'int_steps': 25,
              'int_times': 10,
              'time_step': 0.01,
              'tau': 0.01,
              'agrid': 0.5,
              'box_l': 9.0,
              'dens': 0.85,
              'viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5,
              'skin': 0.2,
              'temp_confidence': 10}
    if espressomd.has_features("SHANCHEN"):
        params.update({"dens": 2 * [params["dens"]]})

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
        self.n_col_part = 1000
        self.system.part.add(pos=np.random.random(
            (self.n_col_part, 3)) * self.params["box_l"])
        if espressomd.has_features("MASS"):
            self.system.part[:].mass = 0.1 + np.random.random(
                len(self.system.part))

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

        # Cache the lb nodes
        lb_nodes = []
        n_nodes = int(self.params['box_l'] / self.params['agrid'])
        for i in range(n_nodes):
            for j in range(n_nodes):
                for k in range(n_nodes):
                    lb_nodes.append(self.lbf[i, j, k])

        # Integration
        for i in range(self.params['int_times']):
            self.system.integrator.run(self.params['int_steps'])
            
            # Summation vars
            fluid_mass = 0.0
            fluid_temp = 0.0
            
            # Go over lb lattice
            for lb_node in lb_nodes:
                fluid_mass += lb_node.density[0]
                fluid_temp += np.sum(lb_node.velocity**2) *lb_node.density[0] 
            
            # Normalize
            fluid_mass/=len(lb_nodes)
            fluid_temp *=  self.system.volume() /(3. * len(lb_nodes)**2)

            # check mass conversation
            self.assertAlmostEqual(fluid_mass,self.params["dens"],delta=self.params["mass_prec_per_node"])

            # check momentum conservation
            np.testing.assert_allclose(
                self.system.analysis.analyze_linear_momentum(),self.tot_mom,
                atol=self.params['mom_prec'])

            # Calc particle temperature
            e = self.system.analysis.energy()
            temp_particle = 2.0 / self.dof * e["kinetic"] / self.n_col_part
           
            # Update lists
            all_temp_particle.append(temp_particle)
            all_temp_fluid.append(fluid_temp)

        #import scipy.stats
        #temp_prec_particle = scipy.stats.norm.interval(0.95, loc=self.params["temp"], scale=np.std(all_temp_particle,ddof=1))[1] -self.params["temp"]
        #temp_prec_fluid = scipy.stats.norm.interval(0.95, loc=self.params["temp"], scale=np.std(all_temp_fluid,ddof=1))[1] -self.params["temp"]
        temp_prec_particle=0.05*self.params["temp"]
        temp_prec_fluid=0.05*self.params["temp"]


        self.assertAlmostEqual(np.mean(all_temp_fluid),self.params["temp"],delta=temp_prec_fluid)
        self.assertAlmostEqual(np.mean(all_temp_particle),self.params["temp"],delta=temp_prec_particle)

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

    @ut.skipIf(not espressomd.has_features("EXTERNAL_FORCES"),
               "Features not available, skipping test!")
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
            np.copy(self.system.part[0].f), -self.params['friction'] * (v_part - v_fluid), atol=1E-6)

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
                np.copy(self.lbf[n].velocity), fluid_velocity, atol=1E-6)


@ut.skipIf(not espressomd.has_features(["LB"]),
           "Features not available, skipping test!")
class TestLBCPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid
        self.params.update({"mom_prec": 1E-9, "mass_prec_per_node": 5E-8})


@ut.skipIf(
    not espressomd.has_features(
        ["LB_GPU"]) or espressomd.has_features('SHANCHEN'),
           "Features not available, skipping test!")
class TestLBGPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU
        self.params.update({"mom_prec": 1E-3, "mass_prec_per_node": 1E-5})


if __name__ == "__main__":
    ut.main()
