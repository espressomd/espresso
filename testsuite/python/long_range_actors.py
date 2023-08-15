#
# Copyright (C) 2021-2022 The ESPResSo project
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

import espressomd
import espressomd.magnetostatics
import espressomd.electrostatics
import espressomd.electrostatic_extensions

import numpy as np
import itertools
import unittest as ut
import unittest_decorators as utx


class Test(ut.TestCase):
    system = espressomd.System(box_l=[20., 20., 20.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    n_nodes = system.cell_system.get_state()["n_nodes"]
    original_node_grid = tuple(system.cell_system.node_grid)

    def setUp(self):
        self.system.box_l = [20., 20., 20.]
        self.system.cell_system.node_grid = self.original_node_grid

    def tearDown(self):
        self.system.part.clear()
        if espressomd.has_features(["ELECTROSTATICS"]):
            self.system.electrostatics.clear()
        if espressomd.has_features(["DIPOLES"]):
            self.system.magnetostatics.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()
        self.system.periodicity = [True, True, True]
        self.system.cell_system.set_regular_decomposition()

    def add_charged_particles(self):
        self.system.part.add(pos=[[0., 0., 0.], [0.5, 0.5, 0.5]], q=[-1., 1.])

    def add_magnetic_particles(self):
        self.system.part.add(pos=[[0.01, 0.01, 0.01], [0.5, 0.5, 0.5]],
                             dip=[(1., 0., 0.), (-1., 0., 0.)],
                             rotation=2 * [(True, True, True)])

    def add_icc_particles(self):
        pos = [[0., 0., 0.], [1., 0., 0.]]
        q = [0.0001, -0.0001]
        p_slice = self.system.part.add(pos=pos, q=q)
        areas = self.system.box_l[0] * self.system.box_l[1] / 2. * np.ones(2)
        normals = 2 * [(0., 0., 1.)]
        return p_slice, normals, areas

    def setup_icc_particles_and_solver(self):
        part_slice, normals, areas = self.add_icc_particles()
        icc = espressomd.electrostatic_extensions.ICC(
            n_icc=len(part_slice),
            normals=normals,
            areas=areas,
            epsilons=np.ones_like(areas),
            first_id=part_slice.id[0],
            max_iterations=100,
            check_neutrality=False,
        )
        return icc, part_slice

    def valid_p3m_parameters(self):
        return {"prefactor": 1., "mesh": 32, "cao": 7, "accuracy": 1e-5,
                "r_cut": 1.25625, "alpha": 1.50505, "tune": False,
                "check_neutrality": False}

    def valid_dp3m_parameters(self):
        return {"prefactor": 1., "mesh": 26, "cao": 7, "accuracy": 1e-4,
                "r_cut": 4.80, "alpha": 0.62626, "tune": False}

    @utx.skipIfMissingFeatures(["P3M"])
    def test_electrostatics_registration(self):
        import espressomd.highlander
        icc, _ = self.setup_icc_particles_and_solver()
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())
        p3m_new = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())

        self.assertIsNone(p3m.call_method("unknown"))
        self.assertIsNone(icc.call_method("unknown"))

        with self.assertRaisesRegex(RuntimeError, "An electrostatics solver is needed by ICC"):
            self.system.electrostatics.extension = icc

        self.system.electrostatics.solver = p3m
        self.system.electrostatics.extension = icc
        if espressomd.has_features(["NPT"]):
            with self.assertRaisesRegex(Exception, "ERROR: ICC does not work in the NPT ensemble"):
                self.system.integrator.set_isotropic_npt(
                    ext_pressure=2., piston=0.01)
                self.system.integrator.run(0)
            self.system.integrator.set_vv()
        with self.assertRaisesRegex(RuntimeError, "Cannot change solver when an extension is active"):
            self.system.electrostatics.solver = p3m_new
        with self.assertRaisesRegex(RuntimeError, "Cannot change solver when an extension is active"):
            self.system.electrostatics.solver = None
        self.assertIsNotNone(self.system.electrostatics.extension)
        self.system.electrostatics.clear()
        self.assertIsNone(self.system.electrostatics.solver)
        self.assertIsNone(self.system.electrostatics.extension)
        self.system.integrator.run(0)
        self.assertIsNone(self.system.electrostatics.call_method("unknown"))

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_magnetostatics_registration(self):
        import espressomd.highlander
        dp3m = espressomd.magnetostatics.DipolarP3M(
            **self.valid_dp3m_parameters())

        self.assertIsNone(dp3m.call_method("unknown"))
        self.system.magnetostatics.solver = dp3m
        self.assertIsNone(self.system.magnetostatics.call_method("unknown"))

    def check_obs_stats(self, key):
        # check observable statistics have the correct shape
        pressure_tensor = self.system.analysis.pressure_tensor()
        pressure_scalar = self.system.analysis.pressure()
        energies_scalar = self.system.analysis.energy()
        for key in [key, (key, 0), (key, 1)]:
            self.assertEqual(np.shape(energies_scalar[key]), ())
            self.assertEqual(np.shape(pressure_scalar[key]), ())
            self.assertEqual(np.shape(pressure_tensor[key]), (3, 3))
            self.assertAlmostEqual(pressure_scalar[key],
                                   np.trace(pressure_tensor[key]) / 3.,
                                   delta=1e-12)
        return pressure_tensor, pressure_scalar

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m_cpu_pressure(self):
        dp3m = espressomd.magnetostatics.DipolarP3M(
            **self.valid_dp3m_parameters())
        self.system.magnetostatics.solver = dp3m
        pressure_tensor, pressure_scalar = self.check_obs_stats("dipolar")
        # DP3M doesn't contribute to the pressure
        np.testing.assert_allclose(pressure_tensor["dipolar"], 0., atol=1e-12)
        np.testing.assert_allclose(pressure_scalar["dipolar"], 0., atol=1e-12)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_cpu_pressure(self):
        self.add_charged_particles()
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())
        self.system.electrostatics.solver = p3m
        self.check_obs_stats("coulomb")

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_gpu_pressure(self):
        self.add_charged_particles()
        p3m = espressomd.electrostatics.P3MGPU(**self.valid_p3m_parameters())
        self.system.electrostatics.solver = p3m
        self.check_obs_stats("coulomb")

    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc_cpu_pressure(self):
        self.add_charged_particles()
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())
        elc = espressomd.electrostatics.ELC(
            actor=p3m, gap_size=2., maxPWerror=1e-3, check_neutrality=False)
        self.system.electrostatics.solver = elc
        pressure_tensor, pressure_scalar = self.check_obs_stats("coulomb")
        # ELC doesn't contribute to the pressure
        pressure_tensor_far_field = pressure_tensor[("coulomb", 1)]
        pressure_scalar_far_field = pressure_scalar[("coulomb", 1)]
        pressure_tensor_near_field = pressure_tensor[("coulomb", 0)]
        pressure_scalar_near_field = pressure_scalar[("coulomb", 0)]
        np.testing.assert_allclose(pressure_tensor_far_field, 0., atol=1e-12)
        np.testing.assert_allclose(pressure_scalar_far_field, 0., atol=1e-12)
        np.testing.assert_allclose(pressure_tensor_near_field, 0., atol=1e-12)
        np.testing.assert_allclose(pressure_scalar_near_field, 0., atol=1e-12)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_rf_pressure(self):
        self.add_charged_particles()
        actor = espressomd.electrostatics.ReactionField(
            prefactor=1., kappa=2., epsilon1=1., epsilon2=2., r_cut=2.)
        self.system.electrostatics.solver = actor
        pressure_tensor, pressure_scalar = self.check_obs_stats("coulomb")
        # actor doesn't contribute to the pressure in the far field
        pressure_tensor_far_field = pressure_tensor[("coulomb", 1)]
        pressure_scalar_far_field = pressure_scalar[("coulomb", 1)]
        np.testing.assert_allclose(pressure_tensor_far_field, 0., atol=1e-12)
        np.testing.assert_allclose(pressure_scalar_far_field, 0., atol=1e-12)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_dh_pressure(self):
        self.add_charged_particles()
        actor = espressomd.electrostatics.DH(prefactor=1., kappa=1., r_cut=1.)
        self.system.electrostatics.solver = actor
        pressure_tensor, pressure_scalar = self.check_obs_stats("coulomb")
        # actor doesn't contribute to the pressure in the far field
        pressure_tensor_far_field = pressure_tensor[("coulomb", 1)]
        pressure_scalar_far_field = pressure_scalar[("coulomb", 1)]
        np.testing.assert_allclose(pressure_tensor_far_field, 0., atol=1e-12)
        np.testing.assert_allclose(pressure_scalar_far_field, 0., atol=1e-12)

    @utx.skipIfMissingFeatures(["DIPOLES"])
    @ut.skipIf(n_nodes > 1, "only runs for 1 MPI rank")
    def test_mdds_cpu_no_magnetic_particles(self):
        self.system.part.add(pos=2 * [[1., 1., 1.]], dip=2 * [[0., 0., 0.]])
        mdds = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=2.)
        self.system.magnetostatics.solver = mdds
        energy = self.system.analysis.energy()
        self.assertAlmostEqual(energy["dipolar"], 0., delta=1e-12)

    def check_p3m_pre_conditions(self, container, class_p3m):
        params = {"prefactor": 1., "accuracy": 1., "r_cut": 1., "alpha": 1.}

        # P3M pre-condition: cao / mesh[i] < 1
        with self.assertRaisesRegex(RuntimeError, "k-space cutoff .+ is larger than half of box dimension"):
            container.solver = class_p3m(cao=6, mesh=6, tune=False, **params)

        # P3M pre-condition: cao / mesh[i] < 2 / n_nodes[i]
        with self.assertRaisesRegex(RuntimeError, "k-space cutoff .+ is larger than local box dimension"):
            self.system.cell_system.node_grid = [self.n_nodes, 1, 1]
            container.solver = class_p3m(cao=7, mesh=8, tune=False, **params)

    @utx.skipIfMissingFeatures(["P3M"])
    @ut.skipIf(n_nodes < 3, "only runs for 3+ MPI ranks")
    def test_p3m_cpu_pre_condition_exceptions(self):
        self.add_charged_particles()
        self.check_p3m_pre_conditions(
            self.system.electrostatics,
            espressomd.electrostatics.P3M)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    @ut.skipIf(n_nodes < 3, "only runs for 3+ MPI ranks")
    def test_p3m_gpu_pre_condition_exceptions(self):
        self.add_charged_particles()
        self.check_p3m_pre_conditions(
            self.system.electrostatics,
            espressomd.electrostatics.P3MGPU)

    @utx.skipIfMissingFeatures(["DP3M"])
    @ut.skipIf(n_nodes < 3, "only runs for 3+ MPI ranks")
    def test_dp3m_cpu_pre_condition_exceptions(self):
        self.add_magnetic_particles()
        self.check_p3m_pre_conditions(
            self.system.magnetostatics,
            espressomd.magnetostatics.DipolarP3M)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_cpu_tuning_accuracy_exception(self):
        self.add_charged_particles()
        p3m = espressomd.electrostatics.P3M(prefactor=1., mesh=8, alpha=60.,
                                            r_cut=1.25625, accuracy=1e-5)
        with self.assertRaisesRegex(RuntimeError, "failed to reach requested accuracy"):
            self.system.electrostatics.solver = p3m
        self.assertFalse(p3m.is_tuned)
        self.assertIsNone(self.system.electrostatics.solver)

    def check_p3m_tuning_errors(self, p3m):
        # set an incompatible combination of thermostat and integrators
        self.system.integrator.set_isotropic_npt(ext_pressure=2., piston=0.01)
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        with self.assertRaisesRegex(RuntimeError, r"tuning failed: an exception was thrown while benchmarking the integration loop"):
            self.system.electrostatics.solver = p3m
        self.assertFalse(p3m.is_tuned)
        self.assertIsNone(self.system.electrostatics.solver)

    @utx.skipIfMissingFeatures(["P3M", "NPT"])
    def test_p3m_cpu_tuning_errors(self):
        self.add_charged_particles()
        p3m = espressomd.electrostatics.P3M(prefactor=1., accuracy=1e-3)
        self.check_p3m_tuning_errors(p3m)

    @utx.skipIfMissingFeatures(["DP3M", "NPT"])
    def test_dp3m_cpu_tuning_errors(self):
        self.add_magnetic_particles()
        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=1., accuracy=1e-3)
        self.check_p3m_tuning_errors(dp3m)

    def check_mmm1d_exceptions(self, mmm1d_class):
        mmm1d = mmm1d_class(prefactor=1., maxPWerror=1e-2)

        # check cell system exceptions
        with self.assertRaisesRegex(Exception, "MMM1D requires the N-square cellsystem"):
            self.system.cell_system.set_regular_decomposition()
            self.system.electrostatics.solver = mmm1d
        self.assertIsNone(self.system.electrostatics.solver)
        self.assertFalse(mmm1d.is_tuned)
        self.system.cell_system.set_n_square()

        # check periodicity exceptions
        for periodicity in itertools.product((True, False), repeat=3):
            if periodicity == (False, False, True):
                continue
            self.system.periodicity = periodicity
            with self.assertRaisesRegex(Exception, r"MMM1D requires periodicity \(False, False, True\)"):
                mmm1d = mmm1d_class(prefactor=1., maxPWerror=1e-2)
                self.system.electrostatics.solver = mmm1d
            self.assertIsNone(self.system.electrostatics.solver)
            self.assertFalse(mmm1d.is_tuned)
        self.system.periodicity = (False, False, True)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_mmm1d_cpu_exceptions(self):
        self.system.periodicity = (False, False, True)
        self.check_mmm1d_exceptions(espressomd.electrostatics.MMM1D)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["MMM1D_GPU"])
    @ut.skipIf(n_nodes > 3, "only runs for 3 or less MPI ranks")
    def test_mmm1d_gpu_exceptions(self):
        # VRAM peak memory usage: 700 MiB on 4 MPI cores, 500 on 3 MPI cores
        self.system.periodicity = (False, False, True)
        self.check_mmm1d_exceptions(espressomd.electrostatics.MMM1DGPU)

        with self.assertRaisesRegex(ValueError, "Parameter 'far_switch_radius' must not be larger than box length"):
            espressomd.electrostatics.MMM1DGPU(
                prefactor=1., maxPWerror=1e-2,
                far_switch_radius=2. * self.system.box_l[2])

    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc_tuning_exceptions(self):
        p3m = espressomd.electrostatics.P3M(**self.valid_p3m_parameters())
        elc = espressomd.electrostatics.ELC(
            actor=p3m,
            gap_size=2.,
            maxPWerror=1e-3,
            delta_mid_top=-1.,
            delta_mid_bot=-1.,
            const_pot=True,
            pot_diff=-3,
            check_neutrality=False,
        )
        self.system.part.add(pos=[0., 0., 0.], q=1.)
        with self.assertRaisesRegex(RuntimeError, "ELC does not currently support non-neutral systems"):
            self.system.electrostatics.solver = elc


if __name__ == "__main__":
    ut.main()
