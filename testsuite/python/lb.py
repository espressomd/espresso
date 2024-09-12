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
              'kinematic_viscosity': 3.0}

    system.periodicity = [True, True, True]
    system.time_step = params['tau']
    system.cell_system.skin = 1.0
    if espressomd.gpu_available():
        system.cuda_init_handle.call_method("set_device_id_per_rank")
    interpolation = False
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def setUp(self):
        self.system.box_l = 3 * [6.0]

    def tearDown(self):
        self.system.lb = None
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
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        self.assertTrue(lbf.is_active)
        self.check_properties(lbf)
        self.system.lb = None

        # deactivated actor
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=1)
        self.system.lb = None
        self.assertFalse(lbf.is_active)
        self.check_properties(lbf)

    def check_properties(self, lbf):
        density = self.params["density"]
        agrid = self.params["agrid"]
        tau = self.params["tau"]
        # check LB object
        self.assertAlmostEqual(lbf.tau, tau, delta=self.atol)
        self.assertAlmostEqual(lbf.agrid, agrid, delta=self.atol)
        self.assertAlmostEqual(lbf.kinematic_viscosity, 3., delta=self.atol)
        self.assertAlmostEqual(lbf.density, density, delta=self.atol)
        self.assertAlmostEqual(lbf.kT, 1.0, delta=self.atol)
        self.assertEqual(lbf.seed, 42)
        self.assertEqual(
            lbf.single_precision,
            self.lb_params["single_precision"])
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), [0., 0., 0.], atol=self.atol)
        lbf.kinematic_viscosity = 2.
        self.assertAlmostEqual(lbf.kinematic_viscosity, 2., delta=self.atol)
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
            lbf.single_precision,
            self.lb_params["single_precision"])
        lbf.kinematic_viscosity = 3.
        self.assertAlmostEqual(lbf.kinematic_viscosity, 3., delta=self.atol)
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
        # check node setters update cached velocities (with precision loss)
        lbnode = lbf[0, 0, 0]
        lbnode.last_applied_force = [0., 0., 0.]
        lbnode.velocity = [0., 0., 0.]
        old_pop = np.copy(lbnode.population)
        old_vel = np.copy(lbnode.velocity)
        old_rho = np.copy(lbnode.density)
        lbnode.velocity = [1., 2., 3.]
        new_pop = np.copy(lbnode.population)
        new_vel = np.copy(lbnode.velocity)
        lbnode.population = old_pop
        np.testing.assert_allclose(
            np.copy(lbnode.velocity), old_vel, atol=self.atol * 20.)
        lbnode.population = new_pop
        np.testing.assert_allclose(
            np.copy(lbnode.velocity), new_vel, atol=self.atol * 20.)
        lbnode.population = old_pop
        lbnode.last_applied_force = [0.4, 0.5, 0.6]
        np.testing.assert_allclose(
            np.copy(lbnode.velocity),
            np.copy([0.4, 0.5, 0.6]) / (agrid / tau * old_rho / 2.),
            atol=self.atol * 20.)
        lbnode.last_applied_force = [0., 0., 0.]
        lbnode.population = old_pop
        # check slice setters update cached velocities (with precision loss)
        lbslice = lbf[0:5, 0:5, 0:5]
        lbslice.last_applied_force = [0., 0., 0.]
        lbslice.velocity = [0., 0., 0.]
        old_pop = np.copy(lbslice.population)
        old_vel = np.copy(lbslice.velocity)
        old_rho = np.copy(lbslice.density)
        lbslice.velocity = [1., 2., 3.]
        new_pop = np.copy(lbslice.population)
        new_vel = np.copy(lbslice.velocity)
        lbslice.population = old_pop
        np.testing.assert_allclose(
            np.copy(lbslice.velocity), old_vel, atol=self.atol * 100.)
        lbslice.population = new_pop
        np.testing.assert_allclose(
            np.copy(lbslice.velocity), new_vel, atol=self.atol * 100.)
        lbslice.population = old_pop
        lbslice.last_applied_force = [0.4, 0.5, 0.6]
        vel2force = 2. * tau / agrid
        np.testing.assert_allclose(
            np.copy(lbslice.velocity),
            np.tile(np.copy([0.4, 0.5, 0.6]), (5, 5, 5, 1)) * vel2force /
            np.tile(old_rho.reshape((5, 5, 5, 1)), (3,)),
            atol=self.atol * 20.)
        lbslice.last_applied_force = [0., 0., 0.]
        lbslice.population = old_pop
        # check node boundary conditions
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
        # momentum update: check density conservation when velocity changes,
        # and velocity conservation when density changes
        node = lbf[1, 2, 3]
        density_old = node.density
        density_new = 0.5
        velocity_old = node.velocity
        velocity_new = [0.01, 0.02, 0.03]
        node.velocity = velocity_new
        np.testing.assert_allclose(np.copy(node.density),
                                   np.copy(density_old), atol=self.atol)
        np.testing.assert_allclose(np.copy(node.velocity),
                                   np.copy(velocity_new), atol=self.atol)
        node.density = density_new
        np.testing.assert_allclose(np.copy(node.density),
                                   np.copy(density_new), atol=self.atol)
        np.testing.assert_allclose(np.copy(node.velocity),
                                   np.copy(velocity_new), atol=self.atol)
        node.density = density_old
        node.velocity = velocity_old
        # check slice matches node
        lbslice = lbf[0:5, 0:5, 0:5]
        np.testing.assert_allclose(
            np.copy(lbslice.population)[1, 2, 3, :],
            np.copy(node.population), atol=self.atol)
        np.testing.assert_allclose(
            np.copy(lbslice.velocity)[1, 2, 3, :],
            np.copy(node.velocity), atol=self.atol)
        np.testing.assert_allclose(
            np.copy(lbslice.pressure_tensor)[1, 2, 3, :],
            np.copy(node.pressure_tensor), atol=self.atol)
        np.testing.assert_allclose(
            np.copy(lbslice.pressure_tensor_neq)[1, 2, 3, :],
            np.copy(node.pressure_tensor_neq), atol=self.atol)
        np.testing.assert_allclose(
            np.copy(lbslice.density)[1, 2, 3],
            np.copy(node.density), atol=self.atol)

    def test_raise_if_read_only(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        for key in {'agrid', 'tau', 'density', 'kT', 'single_precision',
                    'shape', 'pressure_tensor', 'seed', 'is_active'}:
            with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) '{key}' is read-only"):
                setattr(lbf, key, 0)

    def test_ctor_exceptions(self):
        def make_kwargs(**kwargs):
            lb_kwargs = {}
            lb_kwargs.update(self.params)
            lb_kwargs.update(self.lb_params)
            lb_kwargs.update(kwargs)
            return lb_kwargs

        with self.assertRaisesRegex(ValueError, "Parameter 'agrid' must be > 0"):
            self.lb_class(**make_kwargs(agrid=0.))
        with self.assertRaisesRegex(ValueError, "Parameter 'agrid' must be > 0"):
            self.lb_class(**make_kwargs(agrid=-1.))
        with self.assertRaisesRegex(ValueError, "Parameter 'tau' must be > 0"):
            self.lb_class(**make_kwargs(tau=0.))
        with self.assertRaisesRegex(ValueError, "Parameter 'density' must be > 0"):
            self.lb_class(**make_kwargs(density=0.))
        with self.assertRaisesRegex(ValueError, "Parameter 'kinematic_viscosity' must be >= 0"):
            self.lb_class(**make_kwargs(kinematic_viscosity=-1.))
        with self.assertRaisesRegex(ValueError, "Parameter 'kT' must be >= 0"):
            self.lb_class(**make_kwargs(kT=-1., seed=42))
        with self.assertRaisesRegex(ValueError, "Parameter 'seed' must be >= 0"):
            self.lb_class(**make_kwargs(kT=0., seed=-42))

    def test_node_exceptions(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        lb_node = lbf[0, 0, 0]
        # check exceptions from LB node
        with self.assertRaisesRegex(RuntimeError, "Property 'boundary_force' is read-only"):
            lb_node.boundary_force = [1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor' is read-only"):
            lb_node.pressure_tensor = np.eye(3, 3)
        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor_neq' is read-only"):
            lb_node.pressure_tensor_neq = np.eye(3, 3)
        with self.assertRaisesRegex(RuntimeError, "Property 'is_boundary' is read-only"):
            lb_node.is_boundary = True
        with self.assertRaisesRegex(NotImplementedError, "Cannot serialize LB fluid node objects"):
            lb_node.__reduce__()
        # check property types
        array_locked = espressomd.utils.array_locked
        self.assertIsInstance(lb_node.pressure_tensor, array_locked)
        self.assertIsInstance(lb_node.pressure_tensor_neq, array_locked)
        # self.assertIsInstance(lb_node.boundary_force, array_locked) # TODO
        self.assertIsInstance(lb_node.velocity, array_locked)
        self.assertIsInstance(lb_node.last_applied_force, array_locked)

    def test_slice_exceptions(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        lb_slice = lbf[:, :, :]
        # check exceptions from LB slice
        with self.assertRaisesRegex(RuntimeError, "Property 'boundary_force' is read-only"):
            lb_slice.boundary_force = [1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor' is read-only"):
            lb_slice.pressure_tensor = np.eye(3, 3)
        with self.assertRaisesRegex(RuntimeError, "Property 'pressure_tensor_neq' is read-only"):
            lb_slice.pressure_tensor_neq = np.eye(3, 3)
        with self.assertRaisesRegex(RuntimeError, "Property 'is_boundary' is read-only"):
            lb_slice.is_boundary = True
        with self.assertRaisesRegex(NotImplementedError, 'Cannot serialize LB fluid slice objects'):
            lb_slice.__reduce__()
        with self.assertRaisesRegex(RuntimeError, "Unknown fluid property 'unknown'"):
            lb_slice.call_method("get_value_shape", name="unknown")
        # check property types
        array_locked = espressomd.utils.array_locked
        self.assertIsInstance(lb_slice.pressure_tensor, array_locked)
        self.assertIsInstance(lb_slice.pressure_tensor_neq, array_locked)
        # self.assertIsInstance(lb_slice.boundary_force, array_locked) # TODO
        self.assertIsInstance(lb_slice.velocity, array_locked)
        self.assertIsInstance(lb_slice.last_applied_force, array_locked)
        # check exceptions from python slices
        with self.assertRaisesRegex(NotImplementedError, "Slices with step != 1 are not supported"):
            lbf[:10:2, :, :]
        with self.assertRaisesRegex(NotImplementedError, "Tuple-based indexing is not supported"):
            lbf[:2, (0, 1), (0, 1)]
        with self.assertRaisesRegex(AttributeError, "Cannot set properties of an empty .+ object"):
            lbf[0:1, 0:1, 0:0].density = []

    def test_lb_slice_set_get(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        ref_density = 1. + np.arange(np.prod(lbf.shape)).reshape(lbf.shape)
        lbf[:, :, :].density = ref_density
        densities = np.copy(lbf[:, :, :].density)
        np.testing.assert_allclose(densities, ref_density, rtol=1e-5)
        self.assertIsNone(lbf[:1, 0, 0].boundary[0])

        # prepare various slicing operations
        slices = []
        for i in range(3):
            slices.append([
                slice(lbf.shape[i]), slice(0, lbf.shape[i]), slice(1, -1),
                slice(0, 0), slice(5, 1), slice(0, -lbf.shape[i] + 1),
                slice(-lbf.shape[i], None), slice(2, lbf.shape[i] - 1), 1])

        # check gettters
        for subset in itertools.product(*slices):
            # skip indexing without any slice
            if not any(isinstance(item, slice) for item in subset):
                continue
            np.testing.assert_allclose(
                np.copy(lbf[subset].density), ref_density[subset], rtol=1e-5)

        # check settters
        for subset in itertools.product(*slices):
            # skip indexing without any slice and skip slices with zero length
            if not any(isinstance(item, slice) for item in subset) or any(
                    isinstance(s, slice) and (s.start or 0) >= (s.stop or 0) for s in subset):
                continue
            lbf[:, :, :].density = ref_density
            lbf[subset].density = -lbf[subset].density
            densities = np.copy(lbf[:, :, :].density)
            np.testing.assert_allclose(
                densities[subset], -ref_density[subset], rtol=1e-5)
            densities[subset] *= -1.
            np.testing.assert_allclose(densities, ref_density, rtol=1e-5)

        # empty slices
        self.assertEqual(lbf[5:2, 0, 0].pressure_tensor.shape, (0, 3, 3))
        self.assertEqual(lbf[5:2, 0, 0].pressure_tensor_neq.shape, (0, 3, 3))
        self.assertEqual(lbf[5:2, 0:0, -1:-1].velocity.shape, (0, 0, 0, 3))

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
        system.lb = lbf
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
        system.lb = None
        with self.assertRaisesRegex(RuntimeError, 'LB not activated'):
            obs.calculate()

    def test_lb_node_set_get(self):
        lbf = self.lb_class(kT=0.0, ext_force_density=[0, 0, 0], **self.params,
                            **self.lb_params)
        self.system.lb = lbf
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
        ext_force_density = np.array([0.1, 0.2, 1.2])
        last_applied_force = np.array([0.2, 0.4, 0.6])
        lbf.ext_force_density = ext_force_density
        node = lbf[1, 2, 3]
        node.velocity = v_fluid
        np.testing.assert_allclose(np.copy(node.velocity), v_fluid, atol=1e-4)
        node.last_applied_force = last_applied_force
        np.testing.assert_allclose(
            np.copy(node.last_applied_force), last_applied_force, atol=1e-4)
        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), ext_force_density, atol=1e-4)

        self.assertEqual(lbf.kT, 0.0)
        self.assertIsNone(lbf.rng_state)
        with self.assertRaisesRegex(RuntimeError, "This LB instance is unthermalized"):
            lbf.rng_state = 5
        with self.assertRaisesRegex(ValueError, "Parameter 'rng_state' must be >= 0"):
            lbf.rng_state = -5

    def test_temperature_mismatch(self):
        self.system.thermostat.set_langevin(kT=2., seed=23, gamma=2.)
        lbf = self.lb_class(kT=1., seed=42, **self.params, **self.lb_params)
        self.system.lb = lbf
        with self.assertRaisesRegex(RuntimeError, "Cannot set parameter 'kT' to 1.0*: there are currently active thermostats with kT=2.0*"):
            self.system.thermostat.set_lb(LB_fluid=lbf, seed=23, gamma=2.)
        self.assertFalse(self.system.thermostat.lb.is_active)
        self.assertTrue(self.system.thermostat.langevin.is_active)

    def test_parameter_change_without_seed(self):
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=23, gamma=2.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=3.0)
        self.assertAlmostEqual(self.system.thermostat.lb.gamma, 3.)
        with self.assertRaisesRegex(RuntimeError, "Temperature change not supported by LB"):
            self.system.thermostat.turn_off()
            self.system.thermostat.set_langevin(kT=2., seed=23, gamma=2.)
        self.assertFalse(self.system.thermostat.langevin.is_active)
        with self.assertRaisesRegex(ValueError, "must be an integer multiple of the MD time_step"):
            self.system.time_step /= 1.7
        self.system.time_step *= 1.
        if espressomd.has_features("ELECTROSTATICS"):
            self.system.electrostatics.solver = espressomd.electrostatics.DH(
                prefactor=1., kappa=1., r_cut=1.)  # should not fail
            self.system.electrostatics.clear()
        with self.assertRaisesRegex(RuntimeError, "MD cell geometry change not supported by LB"):
            self.system.box_l = [1., 2., 3.]
        np.testing.assert_allclose(np.copy(self.system.box_l), 6., atol=1e-7)
        with self.assertRaisesRegex(RuntimeError, "MPI topology change not supported by LB"):
            self.system.cell_system.node_grid = self.system.cell_system.node_grid

    def test_grid_index(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        # check ranges and out-of-bounds access
        shape = lbf.shape
        int_shape = [int(x) for x in shape]  # cast away numpy integer types
        for i in range(3):
            n = [0, 0, 0]
            n[i] -= shape[i]
            lbf[n[0], n[1], n[2]].velocity
            self.assertEqual(lbf[tuple(n)], lbf[0, 0, 0])
            for offset in (int_shape[i] + 1, -(int_shape[i] + 1)):
                n = [0, 0, 0]
                n[i] += offset
                err_msg = rf"provided index \[{str(n)[1:-1]}\] is out of range for shape \[{str(int_shape)[1:-1]}\]"  # nopep8
                with self.assertRaisesRegex(IndexError, err_msg):
                    lbf[tuple(n)].velocity
        # node index
        node = lbf[1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
            node.index = [2, 4, 6]
        np.testing.assert_array_equal(node.index, [1, 2, 3])
        np.testing.assert_array_equal(
            lbf[-1, -1, -1].index, np.array(shape) - 1)

    def test_incompatible_agrid(self):
        """
        LB lattice initialization must raise an exception when either box_l or
        local_box_l aren't integer multiples of agrid.
        """
        with self.assertRaisesRegex(RuntimeError, "Box length not commensurate with agrid"):
            params = self.params.copy()
            params['agrid'] += 1e-6
            self.lb_class(**params, **self.lb_params)

    def test_agrid_rounding(self):
        """Tests agrid*n ~= box_l for a case where rounding down is needed"""
        system = self.system
        old_l = system.box_l

        n_part = 1000
        phi = 0.05
        lj_sig = 1.0
        l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3 / phi)**(1. / 3.)
        system.box_l = [l] * 3 * np.array(system.cell_system.node_grid)
        lbf = self.lb_class(agrid=l / 31, density=1, kinematic_viscosity=1, kT=0,
                            tau=system.time_step, **self.lb_params)
        system.lb = lbf
        system.integrator.run(steps=1)
        system.lb = None
        system.box_l = old_l

    def test_bool_operations_on_node(self):
        lbf = self.lb_class(kT=1.0, seed=42, **self.params, **self.lb_params)
        self.system.lb = lbf
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
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        particle_rtol = 1e-10
        fluid_rtol = 1e-10
        if self.lb_params["single_precision"]:
            particle_rtol = 4e-6
            fluid_rtol = 4e-6

        # Random velocities
        lbf[:, :, :].velocity = np.random.random((*lbf.shape, 3))
        # Test several particle positions
        for pos in ([0, 0, 0], self.system.box_l, self.system.box_l / 2,
                    self.system.box_l / 2 - self.params['agrid'] / 2):
            p = self.system.part.add(pos=pos, v=[1, 2, 3])
            v_part = np.copy(p.v)

            # In the first time step after a system change, LB coupling forces
            # are ignored. Hence, the coupling position is shifted
            coupling_pos = p.pos + self.system.time_step * p.v
            v_fluid = np.copy(lbf.get_interpolated_velocity(pos=coupling_pos))
            # Nodes to which forces will be interpolated
            lb_nodes = tests_common.get_lb_nodes_around_pos(coupling_pos, lbf)

            self.system.integrator.run(1)
            # Check friction force
            np.testing.assert_allclose(
                np.copy(p.f), -self.gamma * (v_part - v_fluid),
                rtol=particle_rtol, atol=0.)

            # check particle/fluid force balance
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), -np.copy(p.f),
                rtol=fluid_rtol, atol=0.)

            # Check that last_applied_force gets cleared
            p.remove()
            self.system.integrator.run(1)
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), [0, 0, 0])

    def test_viscous_coupling_pairs(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)

        particle_rtol = 1e-10
        fluid_rtol = 1e-10
        if self.lb_params["single_precision"]:
            particle_rtol = 4e-6
            fluid_rtol = 4e-6

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
                np.copy(p1.f), -self.gamma * (v_part1 - v_fluid1),
                rtol=particle_rtol, atol=0.)
            np.testing.assert_allclose(
                np.copy(p2.f), -self.gamma * (v_part2 - v_fluid2),
                rtol=particle_rtol, atol=0.)

            # check particle/fluid force balance
            applied_forces = np.array(
                [n.last_applied_force for n in all_coupling_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), -np.copy(p1.f) - np.copy(p2.f),
                rtol=fluid_rtol, atol=0.)

            # Check that last_applied_force gets cleared
            self.system.part.clear()
            self.system.integrator.run(1)
            applied_forces = np.array(
                [n.last_applied_force for n in all_coupling_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), [0, 0, 0])

    def test_viscous_coupling_rounding(self):
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)
        p = self.system.part.add(pos=[-1E-30] * 3, v=[-1, 0, 0])
        self.system.integrator.run(1)
        for _ in range(20):
            self.system.integrator.run(1)
            self.assertTrue(np.all(p.f != 0.0))

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_INERTIALESS_TRACERS")
    def test_tracers_coupling_rounding(self):
        import espressomd.propagation
        lbf = self.lb_class(**self.params, **self.lb_params)
        self.system.lb = lbf
        ext_f = np.array([0.01, 0.02, 0.03])
        lbf.ext_force_density = ext_f
        self.system.thermostat.set_lb(LB_fluid=lbf, seed=3, gamma=self.gamma)
        rtol = self.rtol
        if lbf.single_precision:
            rtol *= 100.
        mode_tracer = espressomd.propagation.Propagation.TRANS_LB_TRACER
        self.system.time = 0.
        p = self.system.part.add(pos=[-1E-30] * 3, propagation=mode_tracer)
        for _ in range(10):
            self.system.integrator.run(1)
            vel = np.copy(p.v) * self.params["density"]
            vel_ref = (self.system.time + lbf.tau * 0.5) * ext_f
            np.testing.assert_allclose(vel, vel_ref, rtol=rtol, atol=0.)

    @ut.skipIf(n_nodes != 2, "test is designed to run on 2 MPI ranks")
    def test_rng(self):
        system = self.system
        system.lb = self.lb_class(kT=15., **self.params, **self.lb_params)
        system.integrator.run(1)
        diff = system.lb[0, :, :].population - system.lb[6, :, :].population
        # if the RNG uses the local cell index instead of the global cell index,
        # the noise will be identical in all blocks, and the RMS is zero
        rms = np.sqrt(np.mean(np.square(diff)))
        self.assertGreater(rms, 0.01, msg="thermal noise might be correlated!")

    def test_thermalization_force_balance(self):
        system = self.system

        system.part.add(pos=np.random.random((1000, 3)) * system.box_l)
        if espressomd.has_features("MASS"):
            system.part.all().mass = 0.1 + np.random.random(len(system.part))

        lbf = self.lb_class(kT=1.5, seed=4, **self.params, **self.lb_params)
        system.lb = lbf
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

        self.system.lb = lbf
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
        self.system.lb = lbf
        n_time_steps = 1
        self.system.integrator.run(n_time_steps)
        # ext_force_density is a force density, therefore v = ext_force_density
        # / dens * tau * (n_time_steps + 0.5)
        fluid_velocity = np.array(ext_force_density) * self.system.time_step * (
            n_time_steps + 0.5) / self.params['density']
        # Check global linear momentum = density * volume * velocity
        rtol = self.rtol
        if hasattr(lbf, "single_precision") and lbf.single_precision:
            rtol *= 10.
        np.testing.assert_allclose(
            np.copy(self.system.analysis.linear_momentum()),
            fluid_velocity * self.params['density'] * self.system.volume(),
            rtol=rtol)
        # Check node velocities
        for node_velocity in lbf[:, :, :].velocity.reshape((-1, 3)):
            np.testing.assert_allclose(
                np.copy(node_velocity), fluid_velocity, atol=1E-6)

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
            kinematic_viscosity=self.params['kinematic_viscosity'],
            density=self.params['density'],
            agrid=self.params['agrid'])

        def params_with_tau(tau):
            params = base_params.copy()
            params.update(tau=tau)
            return params

        lbf = self.lb_class(**params_with_tau(self.system.time_step),
                            **self.lb_params)
        sim_time = 100 * self.params['tau']
        self.system.lb = lbf
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=0.1)
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        probe_pos = np.array(self.system.box_l) / 2.
        v1 = np.copy(lbf.get_interpolated_velocity(pos=probe_pos))
        f1 = np.copy(p.f)
        self.system.lb = None
        # get fresh LBfluid and change time steps
        with self.assertRaises(Exception):
            self.system.lb = \
                self.lb_class(**params_with_tau(0.5 * self.system.time_step),
                              **self.lb_params)
        self.system.lb = None
        with self.assertRaises(Exception):
            self.system.lb = \
                self.lb_class(**params_with_tau(1.1 * self.system.time_step),
                              **self.lb_params)
        self.system.lb = None

        self.system.lb = \
            self.lb_class(**params_with_tau(self.system.time_step),
                          **self.lb_params)

        with self.assertRaisesRegex(ValueError, r"LB tau \(0\.0100[0-9]+\) must be >= MD time_step \(0\.0200[0-9]+\)"):
            self.system.time_step = 2.0 * lbf.get_params()["tau"]
        with self.assertRaisesRegex(ValueError, r"LB tau \(0\.0100[0-9]+\) must be an integer multiple of the MD time_step \(0\.0080[0-9]+\)"):
            self.system.time_step = 0.8 * lbf.get_params()["tau"]

        self.system.lb = None
        self.system.time_step = 0.5 * self.params['tau']
        lbf = self.lb_class(**params_with_tau(self.system.time_step),
                            **self.lb_params)
        self.system.lb = lbf
        self.system.integrator.run(
            int(round(sim_time / self.system.time_step)))
        v2 = np.copy(lbf.get_interpolated_velocity(pos=probe_pos))
        f2 = np.copy(p.f)
        np.testing.assert_allclose(v1, v2, rtol=1e-2)
        np.testing.assert_allclose(f1, f2, rtol=1e-2)


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberlaDoublePrecisionCPU(LBTest, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures("WALBERLA")
class LBTestWalberlaSinglePrecisionCPU(LBTest, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": True}
    atol = 1e-7
    rtol = 5e-5


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBTestWalberlaDoublePrecisionGPU(LBTest, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBTestWalberlaSinglePrecisionGPU(LBTest, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_lattice_class = espressomd.lb.LatticeWalberla
    lb_params = {"single_precision": True}
    atol = 1e-6
    rtol = 2e-4


if __name__ == "__main__":
    ut.main()
