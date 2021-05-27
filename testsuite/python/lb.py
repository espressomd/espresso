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
from copy import copy

import espressomd
import espressomd.lb
from espressomd.observables import LBFluidPressureTensor
import sys

from tests_common import get_lb_nodes_around_pos


class TestLB:

    """
    Basic tests of the lattice-Boltzmann implementation

    * temperature
    * particle viscous coupling
    * application of external force densities
    * setting and retrieving lb node velocities

    """
    system = espressomd.System(box_l=3 * [6.0])
    np.random.seed(1)
    params = {'time_step': 0.01,
              'tau': 0.01,
              'agrid': 0.5,
              'dens': 0.85,
              'viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5}

    system.periodicity = [1, 1, 1]
    system.time_step = params['time_step']
    system.cell_system.skin = 1.0
    lbf = None
    interpolation = False

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def test_properties(self):
        self.lbf = self.lb_class(
            kT=1.0, seed=42, visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step)
        self.system.actors.add(self.lbf)
        self.lbf.seed = 56
        self.system.integrator.run(1)
        self.assertEqual(self.lbf.seed, 57)
        self.lbf[0, 0, 0].velocity = [1, 2, 3]
        np.testing.assert_allclose(
            np.copy(self.lbf[0, 0, 0].velocity), [1, 2, 3], atol=1E-10)
        with self.assertRaises(Exception):
            self.lbf[0, 0, 0].velocity = [1, 2]
        with self.assertRaises(Exception):
            self.lbf[0, 1].velocity = [1, 2, 3]

    def test_raise_if_not_active(self):
        lbf = self.lb_class(visc=1.0, dens=1.0, agrid=1.0, tau=0.1)
        with self.assertRaises(RuntimeError):
            lbf.seed = 2

    def test_pressure_tensor_observable(self):
        """
        Checks agreement between the LBFluidPressureTensor observable and
        per-node pressure tensor summed up over the entire fluid.

        """
        system = self.system
        self.n_col_part = 1000
        system.part.add(
            pos=np.random.random((self.n_col_part, 3)) * self.system.box_l[0],
            v=np.random.random((self.n_col_part, 3)))
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
        pressure_tensor = np.zeros((3, 3))
        agrid = self.params["agrid"]
        for n in self.lbf.nodes():
            pressure_tensor += n.pressure_tensor

        pressure_tensor /= system.volume() / agrid**3

        obs = LBFluidPressureTensor()
        obs_pressure_tensor = obs.calculate()
        np.testing.assert_allclose(
            pressure_tensor,
            obs_pressure_tensor,
            atol=1E-7)
        np.testing.assert_allclose(
            np.copy(self.lbf.pressure_tensor),
            obs_pressure_tensor,
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
        self.assertAlmostEqual(
            self.lbf[0, 0, 0].density,
            self.params['dens'],
            delta=1e-4)

        shape_ref = np.copy(self.system.box_l) / self.params['agrid']
        np.testing.assert_array_equal(self.lbf.shape, shape_ref.astype(int))

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
        out_of_bounds = int(max(self.system.box_l) / self.params['agrid']) + 1
        with self.assertRaises(ValueError):
            _ = self.lbf[out_of_bounds, 0, 0].velocity
        with self.assertRaises(ValueError):
            _ = self.lbf[0, out_of_bounds, 0].velocity
        with self.assertRaises(ValueError):
            _ = self.lbf[0, 0, out_of_bounds].velocity

    def test_incompatible_agrid(self):
        """
        LB lattice initialization must raise an exception when either box_l or
        local_box_l aren't integer multiples of agrid.
        """
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'] + 1e-6,
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
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        print("box_l", self.system.box_l)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            seed=3,
            gamma=self.params['friction'])

        # Random velocities
        for n in self.lbf.nodes():
            n.velocity = np.random.random(3) - .5
        # Test several particle positions
        for pos in ([0, 0, 0], self.system.box_l, self.system.box_l / 2,
                    self.system.box_l / 2 - self.params['agrid'] / 2):
            p = self.system.part.add(pos=pos, v=[1, 2, 3])

            v_part = p.v 
            # In the first time step after a system change, LB coupling forces
            # are ignored. Hence, the coupling position is shifted 
            coupling_pos = p.pos + self.system.time_step * p.v
            v_fluid = self.lbf.get_interpolated_velocity(coupling_pos)
            # Nodes to which forces will be interpolated
            lb_nodes = get_lb_nodes_around_pos(
                coupling_pos, self.lbf)

            self.system.integrator.run(1)
            # Check friction force
            np.testing.assert_allclose(
                np.copy(p.f), -self.params['friction'] * (v_part - v_fluid), atol=1E-10)

            # check particle/fluid force balance
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), -np.copy(p.f), atol=1E-10)

            # Check that last_applied_force gets cleared
            p.remove()
            self.system.integrator.run(1)
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), [0, 0, 0])

    def test_viscous_coupling_pairs(self):
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            seed=3,
            gamma=self.params['friction'])

        # Random velocities
        for n in self.lbf.nodes():
            n.velocity = np.random.random(3) - .5
        # Test several particle positions
        agrid = self.params['agrid']
        offset = -0.99 * np.array((agrid, agrid, agrid))
        for pos in ([agrid / 2, agrid / 2, agrid / 2], self.system.box_l, self.system.box_l / 2,
                    self.system.box_l / 2 - self.params['agrid'] / 2):
            p1 = self.system.part.add(pos=pos, v=[1, 2, 3])
            p2 = self.system.part.add(pos=pos + offset, v=[-2, 1, 0.3])

            v_part1 = p1.v 
            v_part2 = p2.v
            # In the first time step after a system change, LB coupling forces
            # are ignored. Hence, the coupling position is shifted 
            coupling_pos1 = p1.pos + self.system.time_step * p1.v
            coupling_pos2 = p2.pos + self.system.time_step * p2.v

            v_fluid1 = self.lbf.get_interpolated_velocity(coupling_pos1)
            v_fluid2 = self.lbf.get_interpolated_velocity(coupling_pos2)
            # Nodes to which forces will be interpolated
            lb_nodes1 = get_lb_nodes_around_pos(
                coupling_pos1, self.lbf)
            lb_nodes2 = get_lb_nodes_around_pos(
                coupling_pos2, self.lbf)

            all_coupling_nodes = [self.lbf[index] for index in set(
                [n.index for n in (lb_nodes1 + lb_nodes2)])]
            self.system.integrator.run(1)
            # Check friction force
            np.testing.assert_allclose(
                np.copy(p1.f), -self.params['friction'] * (v_part1 - v_fluid1), atol=1E-10)
            np.testing.assert_allclose(
                np.copy(p2.f), -self.params['friction'] * (v_part2 - v_fluid2), atol=1E-10)

            # check particle/fluid force balance
            applied_forces = np.array(
                [n.last_applied_force for n in all_coupling_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), -np.copy(p1.f) - np.copy(p2.f), atol=1E-10)

            # Check that last_applied_force gets cleared
            self.system.part.clear()
            self.system.integrator.run(1)
            applied_forces = np.array(
                [n.last_applied_force for n in all_coupling_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), [0, 0, 0])

    def test_thermalization_force_balance(self):
        system = self.system

        self.system.part.add(
            pos=np.random.random((1000, 3)) * self.system.box_l)
        if espressomd.has_features("MASS"):
            self.system.part[:].mass = 0.1 + np.random.random(
                len(self.system.part))

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

        for _ in range(20):
            system.integrator.run(1)
            particle_force = np.sum(system.part[:].f, axis=0)
            fluid_force = np.sum(
                np.array([n.last_applied_force for n in self.lbf.nodes()]), axis=0)
            np.testing.assert_allclose(particle_force, -fluid_force)

    def test_force_interpolation(self):
        system = self.system
        self.lbf = self.lb_class(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0])

        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            seed=3,
            gamma=self.params['friction'])
        lattice_speed = self.lbf.agrid / self.lbf.tau

        position = np.array([1., 2., 3.])
        position_lb_units = position / self.lbf.agrid
        force = np.array([4., -5., 6.])
        force_lb_units = force / lattice_speed * system.time_step
        self.lbf.add_force_at_pos(position, force_lb_units)

        system.integrator.run(1)

        # the force should be split equally across the 8 nearest vertices
        n_couplings = 0
        for n in self.lbf.nodes():
            if np.sum(np.abs(n.last_applied_force)):
                fluid_force = np.copy(n.last_applied_force)
                np.testing.assert_allclose(fluid_force, force / 8.)
                distance = np.linalg.norm(n.index - position_lb_units)
                self.assertLessEqual(int(np.round(distance**2)), 3)
                n_couplings += 1
        self.assertEqual(n_couplings, 8)

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
        n_time_steps = 1
        self.system.integrator.run(n_time_steps)
        # ext_force_density is a force density, therefore v = ext_force_density
        # / dens * tau * (n_time_steps + 0.5)
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps + 0.5) / self.params['dens']
        # Chck global linear momentum = density * volume * velocity
        np.testing.assert_allclose(
            self.system.analysis.linear_momentum(),
            fluid_velocity *
            self.params['dens'] *
            self.system.volume())

        # Check node velocities
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
        p = self.system.part.add(pos=[0.1, 0.2, 0.3], fix=[1, 1, 1])
        base_params = {}
        base_params.update(
            ext_force_density=[2.3, 1.2, 0.1],
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'])

        def params_with_tau(tau):
            params = copy(base_params)
            params.update(tau=tau)
            return params

        lbf = self.lb_class(**params_with_tau(self.system.time_step))
        sim_time = 100 * self.params['time_step']
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=0.1)
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        probe_pos = np.array(self.system.box_l) / 2.
        v1 = np.copy(lbf.get_interpolated_velocity(probe_pos))
        f1 = np.copy(p.f)
        self.system.actors.clear()
        # get fresh LBfluid and change time steps
        with self.assertRaises(ValueError):
            self.system.actors.add(
                self.lb_class(**params_with_tau(0.5 * self.system.time_step)))
        with self.assertRaises(ValueError):
            self.system.actors.add(
                self.lb_class(params_with_tau(1.1 * self.system.time_step)))

        self.system.actors.clear()
        self.system.actors.add(
            self.lb_class(**params_with_tau(self.system.time_step)))

        with self.assertRaises(Exception):
            self.system.time_step = 2. * lbf.get_params()["tau"]
            self.system.integrator.run(1)

        with self.assertRaises(Exception):
            self.system.time_step = 0.8 * lbf.get_params()["tau"]
        self.system.actors.clear()
        self.system.time_step = 0.5 * self.params['time_step']
        self.system.actors.add(
            self.lb_class(**params_with_tau(self.system.time_step)))
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        v2 = np.copy(lbf.get_interpolated_velocity(probe_pos))
        f2 = np.copy(p.f)
        np.testing.assert_allclose(v1, v2, rtol=1e-2)
        np.testing.assert_allclose(f1, f2, rtol=1e-2)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class TestLBWalberla(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidWalberla

    def test_stress_tensor(self):
        print("stress tensor not implemented for Walberla. skipping test.")

    def test_pressure_tensor_observable(self):
        print("Not supported by Walberla")


if __name__ == "__main__":
    ut.main()
