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
import itertools

import espressomd
import espressomd.lb
import espressomd.observables
import sys
import tests_common


class LBTest:

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
              'density': 0.85,
              'viscosity': 3.0,
              'friction': 2.0,
              'temp': 1.5,
              'gamma': 1.5}

    system.periodicity = [1, 1, 1]
    system.time_step = params['time_step']
    system.cell_system.skin = 1.0
    interpolation = False

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params['time_step']

    def test_properties(self):
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        # check property getters work before the object is initialized
        self.assertAlmostEqual(lbf.tau, self.system.time_step, delta=self.atol)
        self.assertAlmostEqual(
            lbf.agrid,
            self.params['agrid'],
            delta=self.atol)
        self.assertAlmostEqual(lbf.viscosity, 3., delta=self.atol)
        self.assertAlmostEqual(lbf.density, 0.85, delta=self.atol)
        self.assertAlmostEqual(lbf.kT, 1.0, delta=self.atol)
        self.assertEqual(lbf.seed, 42)
        self.assertEqual(
            lbf.is_single_precision,
            self.lb_params['single_precision'])
        self.assertFalse(lbf.is_active)
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), [0., 0., 0.], atol=self.atol)
        # check property setters work before the object is initialized
        lbf.viscosity = 2.
        self.assertAlmostEqual(lbf.viscosity, 2., delta=self.atol)
        ext_f = [0.01, 0.02, 0.03]
        lbf.ext_force_density = ext_f
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), ext_f, atol=self.atol)
        # check property getters/setters work after the object is initialized
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        self.assertEqual(lbf.rng_state, 0)
        self.system.integrator.run(1)
        self.assertEqual(lbf.rng_state, 1)
        lbf.rng_state = 56
        self.system.integrator.run(1)
        self.assertEqual(lbf.rng_state, 57)
        self.assertAlmostEqual(lbf.tau, self.system.time_step, delta=self.atol)
        self.assertAlmostEqual(
            lbf.agrid,
            self.params['agrid'],
            delta=self.atol)
        self.assertAlmostEqual(lbf.kT, 1.0, delta=self.atol)
        self.assertEqual(lbf.seed, 42)
        self.assertEqual(
            lbf.is_single_precision,
            self.lb_params['single_precision'])
        self.assertTrue(lbf.is_active)
        lbf.viscosity = 3.
        self.assertAlmostEqual(lbf.viscosity, 3., delta=self.atol)
        ext_force_density = [0.02, 0.05, 0.07]
        lbf.ext_force_density = ext_force_density
        np.testing.assert_allclose(np.copy(lbf.ext_force_density),
                                   ext_force_density, atol=self.atol)
        # check node getters/setters
        lbf[0, 0, 0].velocity = [1, 2, 3]
        np.testing.assert_allclose(
            np.copy(lbf[0, 0, 0].velocity), [1, 2, 3], atol=self.atol)
        with self.assertRaises(Exception):
            lbf[0, 0, 0].velocity = [1, 2]
        with self.assertRaises(Exception):
            lbf[0, 1].velocity = [1, 2, 3]

    def test_raise_if_read_only(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        for key in {'agrid', 'tau', 'density', 'kT', 'is_single_precision',
                    'shape', 'pressure_tensor', 'seed', 'is_active'}:
            with self.assertRaisesRegex(RuntimeError, f"Parameter '{key}' is read-only"):
                setattr(lbf, key, 0)

    def test_raise_if_not_active(self):
        lbf = self.lb_class(**self.params, **self.lb_params)

        # check exceptions from LB actor
        with self.assertRaisesRegex(RuntimeError, "Cannot set 'rng_state' before walberla is initialized"):
            lbf.rng_state = 5
        with self.assertRaisesRegex(RuntimeError, "Cannot get 'rng_state' before walberla is initialized"):
            lbf.rng_state
        with self.assertRaisesRegex(RuntimeError, "LB not activated"):
            lbf.pressure_tensor
        with self.assertRaisesRegex(RuntimeError, "LB not activated"):
            lbf.get_interpolated_velocity([0, 0, 0])

        # check exceptions from LB node
        self.system.actors.add(lbf)
        node = lbf[0, 0, 0]
        self.system.actors.remove(lbf)
        with self.assertRaises(RuntimeError):
            node.density
        with self.assertRaises(RuntimeError):
            node.density = 1.
        with self.assertRaises(RuntimeError):
            node.velocity
        with self.assertRaises(RuntimeError):
            node.velocity = [1, 1, 1]
        with self.assertRaises(RuntimeError):
            node.boundary_force
        with self.assertRaises(RuntimeError):
            node.boundary
        with self.assertRaises(RuntimeError):
            node.last_applied_force
        with self.assertRaises(RuntimeError):
            node.last_applied_force = [1, 1, 1]
        with self.assertRaises(RuntimeError):
            node.pressure_tensor
        with self.assertRaises(NotImplementedError):
            node.pressure_tensor = np.eye(3, 3)
        with self.assertRaises(RuntimeError):
            node.is_boundary
        with self.assertRaises(NotImplementedError):
            node.is_boundary = 1
        with self.assertRaises(RuntimeError):
            node.population
        with self.assertRaises(RuntimeError):
            node.population = np.zeros(19)

    def test_pressure_tensor_observable(self):
        """
        Checks agreement between the ``LBFluidPressureTensor`` observable and
        per-node pressure tensor summed up over the entire fluid.

        """
        system = self.system
        self.n_col_part = 1000
        system.part.add(
            pos=np.random.random((self.n_col_part, 3)) * self.system.box_l[0],
            v=np.random.random((self.n_col_part, 3)))
        system.thermostat.turn_off()

        lbf = self.lb_class(kT=1., seed=1, ext_force_density=[0, 0, 0],
                            **self.params, **self.lb_params)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        system.integrator.run(10)
        pressure_tensor = np.zeros((3, 3))
        agrid = self.params["agrid"]
        for n in lbf.nodes():
            pressure_tensor += n.pressure_tensor

        pressure_tensor /= system.volume() / agrid**3

        obs = espressomd.observables.LBFluidPressureTensor()
        obs_pressure_tensor = obs.calculate()
        np.testing.assert_allclose(
            pressure_tensor, obs_pressure_tensor, atol=1E-7)
        np.testing.assert_allclose(
            np.copy(lbf.pressure_tensor), obs_pressure_tensor, atol=1E-10)

    def test_lb_node_set_get(self):
        lbf = self.lb_class(kT=0.0, ext_force_density=[0, 0, 0], **self.params,
                            **self.lb_params)
        self.system.actors.add(lbf)
        self.assertAlmostEqual(
            lbf[0, 0, 0].density, self.params['density'], delta=1e-4)

        shape_ref = np.copy(self.system.box_l) / self.params['agrid']
        np.testing.assert_array_equal(lbf.shape, shape_ref.astype(int))

        v_fluid = np.array([1.2, 4.3, 0.2])
        lbf[0, 0, 0].velocity = v_fluid
        np.testing.assert_allclose(
            np.copy(lbf[0, 0, 0].velocity), v_fluid, atol=1e-4)
        density = 0.234
        lbf[0, 0, 0].density = density
        self.assertAlmostEqual(lbf[0, 0, 0].density, density, delta=1e-4)

        self.assertEqual(lbf[3, 2, 1].index, (3, 2, 1))
        ext_force_density = [0.1, 0.2, 1.2]
        lbf.ext_force_density = ext_force_density
        lbf[1, 2, 3].velocity = v_fluid
        np.testing.assert_allclose(
            np.copy(lbf[1, 2, 3].velocity), v_fluid, atol=1e-4)
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), ext_force_density, atol=1e-4)

        self.assertEqual(lbf.kT, 0.0)
        rng_error_msg = 'The LB does not use a random number generator'
        with self.assertRaisesRegex(RuntimeError, rng_error_msg):
            lbf.rng_state
        with self.assertRaisesRegex(RuntimeError, rng_error_msg):
            lbf.rng_state = 5

    def test_parameter_change_without_seed(self):
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=23, gamma=2.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=3.0)

    def test_grid_index(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.actors.add(lbf)
        out_of_bounds = int(max(self.system.box_l) / self.params['agrid']) + 1
        with self.assertRaises(ValueError):
            lbf[out_of_bounds, 0, 0].velocity
        with self.assertRaises(ValueError):
            lbf[0, out_of_bounds, 0].velocity
        with self.assertRaises(ValueError):
            lbf[0, 0, out_of_bounds].velocity

    def test_incompatible_agrid(self):
        """
        LB lattice initialization must raise an exception when either box_l or
        local_box_l aren't integer multiples of agrid.
        """
        lbf = self.lb_class(
            viscosity=self.params['viscosity'],
            density=self.params['density'],
            agrid=self.params['agrid'] + 1e-6,
            tau=self.system.time_step,
            ext_force_density=[0, 0, 0],
            **self.lb_params)
        print("\nTesting LB error messages:", file=sys.stderr)
        sys.stderr.flush()
        with self.assertRaises(Exception):
            self.system.actors.add(lbf)
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
        lbf = self.lb_class(agrid=l / 31, density=1, viscosity=1, kT=0,
                            tau=system.time_step, **self.lb_params)
        system.actors.add(lbf)
        system.integrator.run(steps=1)
        system.actors.clear()
        system.box_l = old_l

    def test_bool_operations_on_node(self):
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.actors.add(lbf)
        # test __eq()__ where a node is equal to itself and not equal to any
        # other node
        assert lbf[0, 0, 0] == lbf[0, 0, 0]
        x, y, z = range(int(self.system.box_l[0])), range(
            int(self.system.box_l[1])), range(int(self.system.box_l[2]))
        nodes = [lbf[i, j, k] for i, j, k in itertools.product(x, y, z)]
        nodes.remove(lbf[0, 0, 0])
        assert all(lbf[0, 0, 0] != node for node in nodes)
        # test __hash()__ intercept to identify nodes based on index rather
        # than name. set() constructor runs hash()
        subset1, subset2 = nodes[:-10], nodes[-10:]
        assert len(set(subset1 + subset1)) == len(subset1)
        assert len(set(subset1 + subset2)) == len(subset1) + len(subset2)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_viscous_coupling(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(
            LB_fluid=lbf, seed=3, gamma=self.params['friction'])

        # Random velocities
        for n in lbf.nodes():
            n.velocity = np.random.random(3) - .5
        # Test several particle positions
        for pos in ([0, 0, 0], self.system.box_l, self.system.box_l / 2,
                    self.system.box_l / 2 - self.params['agrid'] / 2):
            p = self.system.part.add(pos=pos, v=[1, 2, 3])

            v_part = p.v
            # In the first time step after a system change, LB coupling forces
            # are ignored. Hence, the coupling position is shifted
            coupling_pos = p.pos + self.system.time_step * p.v
            v_fluid = lbf.get_interpolated_velocity(coupling_pos)
            # Nodes to which forces will be interpolated
            lb_nodes = tests_common.get_lb_nodes_around_pos(coupling_pos, lbf)

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
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(
            LB_fluid=lbf, seed=3, gamma=self.params['friction'])

        # Random velocities
        for n in lbf.nodes():
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

            v_fluid1 = lbf.get_interpolated_velocity(coupling_pos1)
            v_fluid2 = lbf.get_interpolated_velocity(coupling_pos2)
            # Nodes to which forces will be interpolated
            lb_nodes1 = tests_common.get_lb_nodes_around_pos(
                coupling_pos1, lbf)
            lb_nodes2 = tests_common.get_lb_nodes_around_pos(
                coupling_pos2, lbf)

            all_coupling_nodes = [lbf[index] for index in set(
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

        system.part.add(pos=np.random.random((1000, 3)) * system.box_l)
        if espressomd.has_features("MASS"):
            system.part[:].mass = 0.1 + np.random.random(len(system.part))

        lbf = self.lb_class(kT=self.params['temp'], seed=4, **self.params,
                            **self.lb_params)
        system.actors.add(lbf)
        system.thermostat.set_lb(
            LB_fluid=lbf, seed=3, gamma=self.params['friction'])

        for _ in range(20):
            system.integrator.run(1)
            particle_force = np.sum(system.part[:].f, axis=0)
            fluid_force = np.sum(
                np.array([n.last_applied_force for n in lbf.nodes()]), axis=0)
            np.testing.assert_allclose(
                particle_force, -fluid_force, rtol=self.rtol)

    def test_force_interpolation(self):
        lbf = self.lb_class(**self.params, **self.lb_params)

        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(
            LB_fluid=lbf, seed=3, gamma=self.params['friction'])
        lattice_speed = lbf.agrid / lbf.tau

        position = np.array([1., 2., 3.])
        position_lb_units = position / lbf.agrid
        force = np.array([4., -5., 6.])
        force_lb_units = force / lattice_speed * self.system.time_step
        lbf.add_force_at_pos(position, force_lb_units)

        self.system.integrator.run(1)

        # the force should be split equally across the 8 nearest vertices
        n_couplings = 0
        for n in lbf.nodes():
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
        lbf = self.lb_class(ext_force_density=ext_force_density, **self.params,
                            **self.lb_params)
        self.system.actors.add(lbf)
        n_time_steps = 1
        self.system.integrator.run(n_time_steps)
        # ext_force_density is a force density, therefore v = ext_force_density
        # / dens * tau * (n_time_steps + 0.5)
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps + 0.5) / self.params['density']
        # Chck global linear momentum = density * volume * velocity
        rtol = self.rtol
        if hasattr(lbf, 'is_single_precision') and lbf.is_single_precision:
            rtol = 2e-4
        np.testing.assert_allclose(
            self.system.analysis.linear_momentum(),
            fluid_velocity * self.params['density'] * self.system.volume(),
            rtol=rtol)

        # Check node velocities
        for n in lbf.nodes():
            np.testing.assert_allclose(
                np.copy(n.velocity), fluid_velocity, atol=1E-6,
                err_msg=f"Fluid node velocity not as expected on node {n.index}")

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
            viscosity=self.params['viscosity'],
            density=self.params['density'],
            agrid=self.params['agrid'])

        def params_with_tau(tau):
            params = base_params.copy()
            params.update(tau=tau)
            return params

        lbf = self.lb_class(**params_with_tau(self.system.time_step),
                            **self.lb_params)
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
        with self.assertRaises(Exception):
            self.system.actors.add(
                self.lb_class(**params_with_tau(0.5 * self.system.time_step),
                              **self.lb_params))
        self.system.actors.clear()
        with self.assertRaises(Exception):
            self.system.actors.add(
                self.lb_class(**params_with_tau(1.1 * self.system.time_step),
                              **self.lb_params))
        self.system.actors.clear()

        self.system.actors.add(
            self.lb_class(**params_with_tau(self.system.time_step),
                          **self.lb_params))

        with self.assertRaises(Exception):
            self.system.time_step = 2. * lbf.get_params()["tau"]
            self.system.integrator.run(1)

        with self.assertRaises(Exception):
            self.system.time_step = 0.8 * lbf.get_params()["tau"]
        self.system.actors.clear()
        self.system.time_step = 0.5 * self.params['time_step']
        self.system.actors.add(
            self.lb_class(**params_with_tau(self.system.time_step),
                          **self.lb_params))
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        v2 = np.copy(lbf.get_interpolated_velocity(probe_pos))
        f2 = np.copy(p.f)
        np.testing.assert_allclose(v1, v2, rtol=1e-2)
        np.testing.assert_allclose(f1, f2, rtol=1e-2)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBTestWalberla(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBTestWalberlaSinglePrecision(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    atol = 1e-7
    rtol = 4e-5


if __name__ == "__main__":
    ut.main()
