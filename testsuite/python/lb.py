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
        self.system.box_l = 3 * [6.0]

    def test_properties(self):
        self.lbf = self.lb_class(
            kT=1.0, seed=42, visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step)
        self.system.actors.add(self.lbf)
        with self.assertRaises(ValueError):
            self.lbf.tau = -0.1
        self.assertAlmostEqual(self.lbf.tau, self.system.time_step, places=6)
        with self.assertRaises(ValueError):
            self.lbf.density = -0.1
        self.lbf.density = 1.0
        with self.assertRaises(ValueError):
            self.lbf.viscosity = -0.1
        self.lbf.density = 2.4
        with self.assertRaises(ValueError):
            self.lbf.density = -2.4
        self.assertAlmostEqual(self.lbf.density, 2.4, places=6)
        self.lbf.seed = 56
        self.system.integrator.run(1)
        self.assertEqual(self.lbf.seed, 57)
        self.lbf.tau = 0.2
        self.assertAlmostEqual(self.lbf.tau, 0.2, places=6)
        with self.assertRaises(ValueError):
            self.lbf.set_params(bulk_visc=-1.2)
        self.lbf.set_params(bulk_visc=1.2)
        self.assertAlmostEqual(
            self.lbf.get_params()['bulk_visc'], 1.2, places=6)
        with self.assertRaises(ValueError):
            self.lbf.set_params(gamma_odd=1.3)
        self.lbf.set_params(gamma_odd=0.3)
        self.assertAlmostEqual(
            self.lbf.get_params()['gamma_odd'], 0.3, places=6)
        with self.assertRaises(ValueError):
            self.lbf.set_params(gamma_even=1.4)
        self.lbf.set_params(gamma_even=0.4)
        self.assertAlmostEqual(
            self.lbf.get_params()['gamma_even'], 0.4, places=6)
        self.lbf[0, 0, 0].velocity = [1, 2, 3]
        np.testing.assert_allclose(
            np.copy(self.lbf[0, 0, 0].velocity), [1, 2, 3], atol=1E-10)
        with self.assertRaises(Exception):
            self.lbf[0, 0, 0].velocity = [1, 2]
        with self.assertRaises(Exception):
            self.lbf[0, 1].velocity = [1, 2, 3]

    def test_raise_if_not_active(self):
        class MockLBFluid(self.lb_class):
            '''LB class mock that ignores runtime errors from agrid and tau.'''
            @property
            def agrid(self):
                return 1.

            @agrid.setter
            def agrid(self, value):
                pass

            @property
            def tau(self):
                return 0.01

            @tau.setter
            def tau(self, value):
                pass

        self.check_raise_if_not_active(self.lb_class, False)
        self.check_raise_if_not_active(MockLBFluid, True)

    def check_raise_if_not_active(self, lb_class, mock):
        lbf = lb_class(visc=1.0, dens=1.0, agrid=1.0, tau=0.1)

        # check exceptions from LB actor
        with self.assertRaises(RuntimeError):
            _ = lbf.density
        with self.assertRaises(RuntimeError):
            lbf.density = 0.2
        with self.assertRaises(RuntimeError):
            _ = lbf.viscosity
        with self.assertRaises(RuntimeError):
            lbf.viscosity = 0.2
        with self.assertRaises(RuntimeError):
            _ = lbf.bulk_viscosity
        with self.assertRaises(RuntimeError):
            lbf.bulk_viscosity = 0.2
        with self.assertRaises(RuntimeError):
            _ = lbf.seed
        with self.assertRaises(RuntimeError):
            lbf.seed = 2
        with self.assertRaises(RuntimeError):
            _ = lbf.kT
        with self.assertRaises(RuntimeError):
            lbf.kT = 2
        with self.assertRaises(RuntimeError):
            _ = lbf.shape
        if not mock:
            with self.assertRaises(RuntimeError):
                _ = lbf.agrid
            with self.assertRaises(RuntimeError):
                lbf.agrid = 0.2
            with self.assertRaises(RuntimeError):
                _ = lbf.tau
            with self.assertRaises(RuntimeError):
                lbf.tau = 0.01
        with self.assertRaises(RuntimeError):
            _ = lbf.pressure_tensor
        with self.assertRaises(NotImplementedError):
            lbf.pressure_tensor = np.eye(3, 3)
        with self.assertRaises(RuntimeError):
            _ = lbf.ext_force_density
        with self.assertRaises(RuntimeError):
            lbf.ext_force_density = [1, 1, 1]
        with self.assertRaises(RuntimeError):
            lbf.get_interpolated_velocity([0, 0, 0])

        # check exceptions from LB node
        self.system.actors.add(lbf)
        node = lbf[0, 0, 0]
        self.system.actors.remove(lbf)
        with self.assertRaises(RuntimeError):
            _ = node.density
        with self.assertRaises(RuntimeError):
            node.density = 1.
        with self.assertRaises(RuntimeError):
            _ = node.velocity
        with self.assertRaises(RuntimeError):
            node.velocity = [1, 1, 1]
        with self.assertRaises(RuntimeError):
            _ = node.pressure_tensor
        with self.assertRaises(NotImplementedError):
            node.pressure_tensor = np.eye(3, 3)
        with self.assertRaises(RuntimeError):
            _ = node.pressure_tensor_neq
        with self.assertRaises(NotImplementedError):
            node.pressure_tensor_neq = np.eye(3, 3)
        with self.assertRaises(RuntimeError):
            _ = node.boundary
        with self.assertRaises(NotImplementedError):
            node.boundary = 1
        with self.assertRaises(RuntimeError):
            _ = node.population
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

        obs = espressomd.observables.LBFluidPressureTensor()
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
        # resize system
        self.system.box_l = self.system.box_l + 1.
        shape_ref = np.copy(self.system.box_l) / self.params['agrid']
        np.testing.assert_array_equal(self.lbf.shape, shape_ref.astype(int))
        np.testing.assert_array_equal(
            np.copy(self.lbf[out_of_bounds, 0, 0].velocity), 0.)

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

    def test_bool_operations_on_node(self):
        self.lbf = self.lb_class(
            kT=1.0, seed=42, visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.system.time_step)
        self.system.actors.add(self.lbf)
        # test __eq()__ where a node is equal to itself and not equal to any
        # other node
        assert self.lbf[0, 0, 0] == self.lbf[0, 0, 0]
        x, y, z = range(int(self.system.box_l[0])), range(
            int(self.system.box_l[1])), range(int(self.system.box_l[2]))
        nodes = [self.lbf[i, j, k] for i, j, k in itertools.product(x, y, z)]
        nodes.remove(self.lbf[0, 0, 0])
        assert all(self.lbf[0, 0, 0] != node for node in nodes)
        # test __hash()__ intercept to identify nodes based on index rather
        # than name. set() constructor runs hash()
        subset1, subset2 = nodes[:-10], nodes[-10:]
        assert len(set(subset1 + subset1)) == len(subset1)
        assert len(set(subset1 + subset2)) == len(subset1) + len(subset2)

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
        p = self.system.part.add(
            pos=[0.5 * self.params['agrid']] * 3, v=v_part, fix=[1, 1, 1])
        self.lbf[0, 0, 0].velocity = v_fluid
        if self.interpolation:
            v_fluid = self.lbf.get_interpolated_velocity(p.pos)
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(p.f), -self.params['friction'] * (v_part - v_fluid), atol=1E-6)

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
        p = self.system.part.add(pos=[0.1, 0.2, 0.3], fix=[1, 1, 1])
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
        f1 = np.copy(p.f)
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
        f2 = np.copy(p.f)
        np.testing.assert_allclose(v1, v2, rtol=1e-5)
        np.testing.assert_allclose(f1, f2, rtol=1e-5)


class TestLBCPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
class TestLBGPU(TestLB, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_viscous_coupling_higher_order_interpolation(self):
        self.interpolation = True
        self.test_viscous_coupling()
        self.interpolation = False


if __name__ == "__main__":
    ut.main()
