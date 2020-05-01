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
from espressomd.observables import LBFluidStress
import sys


class TestLB:

    """
    Basic tests of the lattice-Boltzmann implementation

    * mass and momentum conservation
    * temperature
    * particle viscous coupling
    * application of external force densities
    * setting and retrieving lb node velocities

    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1)
    params = {'int_steps': 15,
              'int_times': 20,
              'time_step': 0.01,
              'tau': 0.01,
              'agrid': 0.5,
              'box_l': 6.0,
              'dens': 0.85,
              'viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5,
              'skin': 1.0,
              'temp_confidence': 10}

    dof = 3.

    system.box_l = [
        params['box_l'],
        params['box_l'],
        params['box_l']]
    system.periodicity = [1, 1, 1]
    system.time_step = params['time_step']
    system.cell_system.skin = params['skin']
    lbf = None
    interpolation = False

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def test_mass_momentum_thermostat(self):
        self.n_col_part = 100
        self.system.part.add(pos=np.random.random(
            (self.n_col_part, 3)) * self.params["box_l"])
        if espressomd.has_features("MASS"):
            self.system.part[:].mass = 0.1 + np.random.random(
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
        for _ in range(self.params['int_times']):
            self.system.integrator.run(self.params['int_steps'])

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

    def test_properties(self):
        self.lbf = self.lb_class(
            kT=1.0, seed=42, visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step)
        self.system.actors.add(self.lbf)
        with self.assertRaises(ValueError):
            self.lbf.density = -0.1
        self.lbf.density = 1.0
        with self.assertRaises(ValueError):
            self.lbf.viscosity = -0.1
        self.density = 2.4
        self.assertEqual(self.density, 2.4)
        self.lbf.seed = 56
        self.system.integrator.run(1)
        self.assertEqual(self.lbf.seed, 57)
        self.lbf.tau = 0.2
        self.assertAlmostEqual(self.lbf.tau, 0.2)

    def test_raise_if_not_active(self):
        lbf = self.lb_class(visc=1.0, dens=1.0, agrid=1.0, tau=0.1)
        with self.assertRaises(RuntimeError):
            lbf.viscosity = 0.2
        with self.assertRaises(RuntimeError):
            lbf.bulk_viscosity = 0.2
        with self.assertRaises(RuntimeError):
            lbf.density = 0.2
        with self.assertRaises(RuntimeError):
            lbf.seed = 2
        with self.assertRaises(RuntimeError):
            lbf.agrid = 0.2

    def test_stress_tensor_observable(self):
        """
        Checks agreement between the LBFluidStress observable and per-node
        stress summed up over the entire fluid.

        """
        system = self.system
        self.n_col_part = 1000
        system.part.add(pos=np.random.random(
            (self.n_col_part, 3)) * self.params["box_l"], v=np.random.random((self.n_col_part, 3)))
        system.thermostat.turn_off()

        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=system.time_step,
            kT=1, ext_force_density=[0, 0, 0], seed=1)
        system.actors.add(self.lbf)
        system.thermostat.set_lb(LB_fluid=self.lbf, seed=1)
        system.integrator.run(10)
        stress = np.zeros((3, 3))
        agrid = self.params["agrid"]
        for n in self.lbf.nodes():
            stress += n.stress

        stress /= system.volume() / agrid**3

        obs = LBFluidStress()
        obs_stress = obs.calculate()
        np.testing.assert_allclose(stress, obs_stress, atol=1E-10)
        np.testing.assert_allclose(
            np.copy(self.lbf.stress),
            obs_stress,
            atol=1E-10)

    def test_lb_node_set_get(self):
        self.lbf = self.lb_class(
            kT=0.0,
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)

        self.assertEqual(self.lbf.shape,
                         (int(self.system.box_l[0] / self.params["agrid"]),
                          int(self.system.box_l[1] / self.params["agrid"]),
                          int(self.system.box_l[2] / self.params["agrid"])))

        v_fluid = np.array([1.2, 4.3, 0.2])
        self.lbf[0, 0, 0].velocity = v_fluid
        np.testing.assert_allclose(
            np.copy(self.lbf[0, 0, 0].velocity), v_fluid, atol=1e-4)
        density = 0.234
        self.lbf[0, 0, 0].density = density
        self.assertAlmostEqual(self.lbf[0, 0, 0].density, density, delta=1e-4)

        self.assertEqual(self.lbf[3, 2, 1].index, (3, 2, 1))
        ext_force_density = [0.1, 0.2, 1.2]
        self.lbf.ext_force_density = ext_force_density
        self.lbf[1, 2, 3].velocity = v_fluid
        np.testing.assert_allclose(
            np.copy(self.lbf[1, 2, 3].velocity),
            v_fluid,
            atol=1e-4)
        np.testing.assert_allclose(
            np.copy(self.lbf.ext_force_density),
            ext_force_density,
            atol=1e-4)

    def test_parameter_change_without_seed(self):
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0],
            kT=1.0,
            seed=42)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=23, gamma=2.0)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=3.0)

    def test_grid_index(self):
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        with self.assertRaises(ValueError):
            _ = self.lbf[
                int(self.params['box_l'] / self.params['agrid']) + 1, 0, 0].velocity
        with self.assertRaises(ValueError):
            _ = self.lbf[
                0, int(self.params['box_l'] / self.params['agrid']) + 1, 0].velocity
        with self.assertRaises(ValueError):
            _ = self.lbf[
                0, 0, int(self.params['box_l'] / self.params['agrid']) + 1].velocity

    def test_incompatible_agrid(self):
        """
        LB lattice initialization must raise an exception when either box_l or
        local_box_l aren't integer multiples of agrid.
        """
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'] + 1e-5,
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        print("\nTesting LB error messages:", file=sys.stderr)
        sys.stderr.flush()
        with self.assertRaises(Exception):
            self.system.actors.add(self.lbf)
        print("End of LB error messages", file=sys.stderr)
        sys.stderr.flush()

    def test_agrid_rounding(self):
        """Tests agrid*n ~= box_l for a case where rounding down is needed"""
        system = self.system
        old_l = system.box_l

        n_part = 1000
        phi = 0.05
        lj_sig = 1.0
        l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3 / phi)**(1. / 3.)
        system.box_l = [l] * 3 * system.cell_system.node_grid
        system.actors.add(self.lb_class(agrid=l / 31, dens=1,
                                        visc=1, kT=0, tau=system.time_step))
        system.integrator.run(steps=1)
        system.actors.clear()
        system.box_l = old_l

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_viscous_coupling(self):
        v_part = np.array([1, 2, 3])
        v_fluid = np.array([1.2, 4.3, 0.2])
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        if self.interpolation:
            self.lbf.set_interpolation_order("quadratic")
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            seed=3,
            gamma=self.params['friction'])
        self.system.part.add(
            pos=[0.5 * self.params['agrid']] * 3, v=v_part, fix=[1, 1, 1])
        self.lbf[0, 0, 0].velocity = v_fluid
        if self.interpolation:
            v_fluid = self.lbf.get_interpolated_velocity(
                self.system.part[0].pos)
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(self.system.part[0].f), -self.params['friction'] * (v_part - v_fluid), atol=1E-6)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_ext_force_density(self):
        ext_force_density = [2.3, 1.2, 0.1]
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=ext_force_density)
        self.system.actors.add(self.lbf)
        n_time_steps = 5
        self.system.integrator.run(n_time_steps)
        # ext_force_density is a force density, therefore v = ext_force_density
        # / dens * tau * (n_time_steps + 0.5)
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps + 0.5) / self.params['dens']
        for n in self.lbf.nodes():
            np.testing.assert_allclose(
                np.copy(n.velocity), fluid_velocity, atol=1E-6, err_msg="Fluid node velocity not as expected on node {}".format(n.index))

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_unequal_time_step(self):
        """
        Checks that LB tau can only be an integer multiple of the MD time_step
        and that different time steps don't affect the physics of a system
        where particles don't move.

        """
        self.system.part.add(pos=[0.1, 0.2, 0.3], fix=[1, 1, 1])
        ext_force_density = [2.3, 1.2, 0.1]
        lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.params['time_step'],
            ext_force_density=ext_force_density,
            kT=0.)
        sim_time = 100 * self.params['time_step']
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=0.1)
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        probe_pos = np.array(self.system.box_l) / 2.
        v1 = np.copy(lbf.get_interpolated_velocity(probe_pos))
        f1 = np.copy(self.system.part[0].f)
        self.system.actors.clear()
        # get fresh LBfluid and change time steps
        lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.params['time_step'],
            ext_force_density=ext_force_density)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=0.1)
        # illegal time_step/ tau combinations
        with self.assertRaises(ValueError):
            lbf.tau = 0.5 * self.system.time_step
        with self.assertRaises(ValueError):
            lbf.tau = 1.1 * self.system.time_step
        with self.assertRaises(ValueError):
            self.system.time_step = 2. * lbf.get_params()["tau"]
        with self.assertRaises(ValueError):
            self.system.time_step = 0.8 * lbf.get_params()["tau"]
        lbf.tau = self.params['time_step']
        self.system.time_step = 0.5 * self.params['time_step']
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        self.system.time_step = self.params['time_step']
        v2 = np.copy(lbf.get_interpolated_velocity(probe_pos))
        f2 = np.copy(self.system.part[0].f)
        np.testing.assert_allclose(v1, v2, rtol=1e-5)
        np.testing.assert_allclose(f1, f2, rtol=1e-5)


class TestLBCPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid
        self.params.update({"mom_prec": 1E-9, "mass_prec_per_node": 5E-8})


@utx.skipIfMissingGPU()
class TestLBGPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU
        self.params.update({"mom_prec": 1E-3, "mass_prec_per_node": 1E-5})

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_viscous_coupling_higher_order_interpolation(self):
        self.interpolation = True
        self.test_viscous_coupling()
        self.interpolation = False


if __name__ == "__main__":
    ut.main()
