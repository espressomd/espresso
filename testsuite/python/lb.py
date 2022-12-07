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
import itertools

import espressomd
import espressomd.lb
import espressomd.utils
import espressomd.observables
import espressomd.electrostatics
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
    gamma = 2.0
    params = {'tau': 0.01,
              'agrid': 0.5,
              'density': 0.85,
              'viscosity': 3.0}

    system.periodicity = [True, True, True]
    system.time_step = params['tau']
    system.cell_system.skin = 1.0
    interpolation = False

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params['tau']

    def test_properties(self):
        # inactive actor
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.assertFalse(lbf.is_active)
        self.check_properties(lbf)

        # activated actor
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        self.assertTrue(lbf.is_active)
        self.check_properties(lbf)
        self.system.actors.remove(lbf)

        # deactivated actor
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        self.system.actors.remove(lbf)
        self.assertFalse(lbf.is_active)
        self.check_properties(lbf)

    def check_properties(self, lbf):
        agrid = self.params['agrid']
        tau = self.system.time_step
        # check LB object
        self.assertAlmostEqual(lbf.tau, tau, delta=self.atol)
        self.assertAlmostEqual(lbf.agrid, agrid, delta=self.atol)
        self.assertAlmostEqual(lbf.viscosity, 3., delta=self.atol)
        self.assertAlmostEqual(lbf.density, 0.85, delta=self.atol)
        self.assertAlmostEqual(lbf.kT, 1.0, delta=self.atol)
        self.assertEqual(lbf.seed, 42)
        self.assertEqual(
            lbf.is_single_precision,
            self.lb_params['single_precision'])
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), [0., 0., 0.], atol=self.atol)
        lbf.viscosity = 2.
        self.assertAlmostEqual(lbf.viscosity, 2., delta=self.atol)
        ext_f = [0.01, 0.02, 0.03]
        lbf.ext_force_density = ext_f
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), ext_f, atol=self.atol)
        self.assertEqual(lbf.rng_state, 0)
        self.system.integrator.run(1)
        self.assertEqual(lbf.rng_state, int(lbf.is_active))
        lbf.rng_state = 56
        self.system.integrator.run(1)
        self.assertEqual(lbf.rng_state, 56 + int(lbf.is_active))
        self.assertAlmostEqual(lbf.tau, tau, delta=self.atol)
        self.assertAlmostEqual(lbf.agrid, agrid, delta=self.atol)
        self.assertAlmostEqual(lbf.kT, 1.0, delta=self.atol)
        self.assertEqual(lbf.seed, 42)
        self.assertEqual(
            lbf.is_single_precision,
            self.lb_params['single_precision'])
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
        with self.assertRaises(RuntimeError):
            lbf[0, 0, 0].velocity = [1, 2]
        with self.assertRaises(TypeError):
            lbf[0, 1].velocity = [1, 2, 3]
        node = lbf[0, 0, 0]
        self.assertIsNone(node.boundary)
        self.assertIsNone(node.boundary_force)
        vbb_ref = espressomd.lb.VelocityBounceBack([1e-6, 2e-6, 3e-6])
        node.boundary = vbb_ref
        np.testing.assert_allclose(
            np.copy(node.boundary.velocity), np.copy(vbb_ref.velocity),
            atol=self.atol)
        with self.assertRaisesRegex(TypeError, "Parameter 'value' must be an instance of VelocityBounceBack or None"):
            node.boundary = vbb_ref.velocity
        # TODO WALBERLA: remove next line (no-op to get code coverage) once
        # the boundary force getter is implemented from the waLBerla side
        self.assertEqual(len(node.boundary_force), 3)

    def test_raise_if_read_only(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        for key in {'agrid', 'tau', 'density', 'kT', 'is_single_precision',
                    'shape', 'pressure_tensor', 'seed', 'is_active'}:
            with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) '{key}' is read-only"):
                setattr(lbf, key, 0)

    def test_ctor_exceptions(self):
        lattice = self.lb_lattice_class(
            agrid=2 * self.params['agrid'], n_ghost_layers=1)
        with self.assertRaisesRegex(ValueError, "cannot provide both 'lattice' and 'agrid'"):
            self.lb_class(lattice=lattice, **self.params, **self.lb_params)
        with self.assertRaisesRegex(ValueError, "density must be a strictly positive number"):
            params = self.params.copy()
            params['density'] = 0.
            self.lb_class(**params, **self.lb_params)
        with self.assertRaisesRegex(ValueError, "kT must be a positive number"):
            params = self.params.copy()
            params['kT'] = -1e-12
            self.lb_class(**params, **self.lb_params)

    def test_node_exceptions(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.actors.add(lbf)
        node = lbf[0, 0, 0]
        # check exceptions from LB node
        with self.assertRaisesRegex(RuntimeError, "Property 'boundary_force' is read-only"):
            node.boundary_force = [1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor' is read-only"):
            node.pressure_tensor = np.eye(3, 3)
        with self.assertRaisesRegex(RuntimeError, "Property 'is_boundary' is read-only"):
            node.is_boundary = True
        with self.assertRaisesRegex(NotImplementedError, 'Cannot serialize LB fluid node objects'):
            node.__reduce__()
        with self.assertRaisesRegex(NotImplementedError, 'Cannot serialize LB fluid slice objects'):
            lbf[0, 0, 0:2].__reduce__()
        # check property types
        array_locked = espressomd.utils.array_locked
        self.assertIsInstance(node.pressure_tensor, array_locked)
        # self.assertIsInstance(node.boundary_force, array_locked) # TODO
        self.assertIsInstance(node.velocity, array_locked)
        self.assertIsInstance(node.last_applied_force, array_locked)

    def test_pressure_tensor_observable(self):
        """
        Checks agreement between the ``LBFluidPressureTensor`` observable and
        per-node pressure tensor summed up over the entire fluid.

        """
        system = self.system
        n_col_part = 1000
        system.part.add(
            pos=np.random.random((n_col_part, 3)) * self.system.box_l[0],
            v=np.random.random((n_col_part, 3)))
        system.thermostat.turn_off()

        lbf = self.lb_class(kT=1., seed=1, ext_force_density=[0, 0, 0],
                            **self.params, **self.lb_params)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        system.integrator.run(10)
        pressure_tensor = np.copy(
            np.mean(lbf[:, :, :].pressure_tensor, axis=(0, 1, 2)))

        obs = espressomd.observables.LBFluidPressureTensor()
        obs_pressure_tensor = obs.calculate()
        np.testing.assert_allclose(
            pressure_tensor, obs_pressure_tensor,
            atol=self.atol, rtol=self.rtol)
        np.testing.assert_allclose(
            np.copy(lbf.pressure_tensor), obs_pressure_tensor,
            atol=1e-12, rtol=1e-12)

        self.assertIsInstance(
            lbf[0, 0, 0].pressure_tensor, espressomd.utils.array_locked)
        self.assertIsInstance(
            lbf.pressure_tensor,
            espressomd.utils.array_locked)
        system.actors.remove(lbf)
        with self.assertRaisesRegex(RuntimeError, 'LB not activated'):
            obs.calculate()

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
        actor = espressomd.electrostatics.DH(prefactor=1., kappa=1., r_cut=1.)
        with self.assertRaisesRegex(Exception, "Temperature change not supported by LB"):
            self.system.thermostat.turn_off()
        with self.assertRaisesRegex(Exception, "Time step change not supported by LB"):
            self.system.time_step /= 2.
        with self.assertRaisesRegex(RuntimeError, "LB does not currently support handling changes of the MD cell geometry"):
            self.system.actors.add(actor)
        self.assertEqual(len(self.system.actors), 1)

    def test_grid_index(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.actors.add(lbf)
        # check ranges and out-of-bounds access
        shape = lbf.shape
        for i in range(3):
            n = [0, 0, 0]
            n[i] -= shape[i]
            lbf[n[0], n[1], n[2]].velocity
            self.assertEqual(lbf[tuple(n)], lbf[0, 0, 0])
            for offset in (shape[i] + 1, -(shape[i] + 1)):
                n = [0, 0, 0]
                n[i] += offset
                err_msg = rf"provided index \[{str(n)[1:-1]}\] is out of range for shape \[{str(list(shape))[1:-1]}\]"
                with self.assertRaisesRegex(IndexError, err_msg):
                    lbf[tuple(n)].velocity
        # node index
        node = lbf[1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
            node.index = [2, 4, 6]
        np.testing.assert_array_equal(node.index, [1, 2, 3])
        retval = node.call_method("override_index", index=[2, 4, 6])
        self.assertEqual(retval, 0)
        np.testing.assert_array_equal(node.index, [2, 4, 6])
        retval = node.call_method("override_index", index=[0, 0, shape[2]])
        self.assertEqual(retval, 1)
        np.testing.assert_array_equal(node.index, [2, 4, 6])
        np.testing.assert_array_equal(
            lbf[-1, -1, -1].index, np.array(shape) - 1)

    def test_incompatible_agrid(self):
        """
        LB lattice initialization must raise an exception when either box_l or
        local_box_l aren't integer multiples of agrid.
        """
        print("\nTesting LB runtime error messages:", file=sys.stderr)
        sys.stderr.flush()
        with self.assertRaises(Exception):
            params = self.params.copy()
            params['agrid'] += 1e-6
            self.lb_class(**params, **self.lb_params)
        print("End of LB runtime error messages", file=sys.stderr)
        sys.stderr.flush()

    def test_agrid_rounding(self):
        """Tests agrid*n ~= box_l for a case where rounding down is needed"""
        system = self.system
        old_l = system.box_l

        n_part = 1000
        phi = 0.05
        lj_sig = 1.0
        l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3 / phi)**(1. / 3.)
        system.box_l = [l] * 3 * np.array(system.cell_system.node_grid)
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
        shape = np.around(self.system.box_l / self.params["agrid"]).astype(int)
        nodes = [
            lbf[ijk] for ijk in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2]))]
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
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        # Random velocities
        lbf[:, :, :].velocity = np.random.random((*lbf.shape, 3))
        # Test several particle positions
        for pos in ([0, 0, 0], self.system.box_l, self.system.box_l / 2,
                    self.system.box_l / 2 - self.params['agrid'] / 2):
            p = self.system.part.add(pos=pos, v=[1, 2, 3])

            v_part = p.v
            # In the first time step after a system change, LB coupling forces
            # are ignored. Hence, the coupling position is shifted
            coupling_pos = p.pos + self.system.time_step * p.v
            v_fluid = lbf.get_interpolated_velocity(pos=coupling_pos)
            # Nodes to which forces will be interpolated
            lb_nodes = tests_common.get_lb_nodes_around_pos(coupling_pos, lbf)

            self.system.integrator.run(1)
            # Check friction force
            np.testing.assert_allclose(
                np.copy(p.f), -self.gamma * (v_part - v_fluid), atol=1E-10)

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
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        # Random velocities
        lbf[:, :, :].velocity = np.random.random((*lbf.shape, 3))
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

            v_fluid1 = lbf.get_interpolated_velocity(pos=coupling_pos1)
            v_fluid2 = lbf.get_interpolated_velocity(pos=coupling_pos2)
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
                np.copy(p1.f), -self.gamma * (v_part1 - v_fluid1), atol=1E-10)
            np.testing.assert_allclose(
                np.copy(p2.f), -self.gamma * (v_part2 - v_fluid2), atol=1E-10)

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
            system.part.all().mass = 0.1 + np.random.random(len(system.part))

        lbf = self.lb_class(kT=1.5, seed=4, **self.params, **self.lb_params)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        for _ in range(20):
            system.integrator.run(1)
            particle_force = np.sum(system.part.all().f, axis=0)
            fluid_force = np.copy(
                np.sum(lbf[:, :, :].last_applied_force, axis=(0, 1, 2)))
            np.testing.assert_allclose(
                particle_force, -fluid_force, rtol=self.rtol)

    def test_force_interpolation(self):
        lbf = self.lb_class(**self.params, **self.lb_params)

        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        position = np.array([1., 2., 3.])
        position_lb_units = position / lbf.agrid
        force = np.array([4., -5., 6.])
        lbf.add_force_at_pos(pos=position, force=force)

        self.system.integrator.run(1)

        # the force should be split equally across the 8 nearest vertices
        n_couplings = 0
        for n in lbf[:, :, :]:
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
        # Walberla TODO: 0.5 needs to be added to time step once f/2
        # correction is back
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps + 0.5) / self.params['density']
        # Chck global linear momentum = density * volume * velocity
        rtol = self.rtol
        ratio = self.system.analysis.linear_momentum() \
            / (fluid_velocity * self.params['density'] * self.system.volume())
        print(ratio, 1 / ratio)
        if hasattr(lbf, 'is_single_precision') and lbf.is_single_precision:
            rtol *= 10.
        # Walberla todo: The factor 1.5 can go, once f/2 correction is back
        np.testing.assert_allclose(
            np.copy(self.system.analysis.linear_momentum()),
            fluid_velocity * self.params['density'] * self.system.volume(),
            rtol=rtol)

        # Check node velocities
        for node_velocity in lbf[:, :, :].velocity.reshape((-1, 3)):
            np.testing.assert_allclose(
                node_velocity, fluid_velocity, atol=1E-6)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_unequal_time_step(self):
        """
        Checks that LB tau can only be an integer multiple of the MD time_step
        and that different time steps don't affect the physics of a system
        where particles don't move.

        """
        p = self.system.part.add(pos=[0.1, 0.2, 0.3], fix=[True, True, True])
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
        sim_time = 100 * self.params['tau']
        self.system.actors.add(lbf)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=0.1)
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        probe_pos = np.array(self.system.box_l) / 2.
        v1 = np.copy(lbf.get_interpolated_velocity(pos=probe_pos))
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

        with self.assertRaisesRegex(ValueError, r"LB tau \(0\.0100[0-9]+\) must be >= MD time_step \(0\.0200[0-9]+\)"):
            self.system.time_step = 2.0 * lbf.get_params()["tau"]
        with self.assertRaisesRegex(ValueError, r"LB tau \(0\.0100[0-9]+\) must be an integer multiple of the MD time_step \(0\.0080[0-9]+\)"):
            self.system.time_step = 0.8 * lbf.get_params()["tau"]

        self.system.actors.clear()
        self.system.time_step = 0.5 * self.params['tau']
        lbf = self.lb_class(**params_with_tau(self.system.time_step),
                            **self.lb_params)
        self.system.actors.add(lbf)
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        v2 = np.copy(lbf.get_interpolated_velocity(pos=probe_pos))
        f2 = np.copy(p.f)
        np.testing.assert_allclose(v1, v2, rtol=1e-2)
        np.testing.assert_allclose(f1, f2, rtol=1e-2)


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberla(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {'single_precision': False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberlaSinglePrecision(LBTest, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {'single_precision': True}
    atol = 1e-7
    rtol = 5e-5


if __name__ == "__main__":
    ut.main()
