#
# Copyright (C) 2013-2022 The ESPResSo project
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
import tests_common

import espressomd.electrostatics


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class Test(ut.TestCase):
    system = espressomd.System(box_l=[10., 10., 10.])
    original_node_grid = system.cell_system.node_grid

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.periodicity = [True, True, True]
        self.system.cell_system.set_regular_decomposition()
        self.system.cell_system.node_grid = self.original_node_grid
        self.system.time_step = 0.01
        self.system.part.add(pos=(0.0, 0.0, 0.0), q=+1.)
        self.system.part.add(pos=(0.1, 0.1, 0.1), q=-1.)

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    if espressomd.has_features(["ELECTROSTATICS"]):
        test_dh = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.DH,
            dict(prefactor=2., kappa=3., r_cut=1.5,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12))

    if espressomd.has_features(["ELECTROSTATICS"]):
        test_rf = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.ReactionField,
            dict(prefactor=2., kappa=3., epsilon1=4., epsilon2=5., r_cut=1.5,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12))

    if espressomd.has_features(["P3M"]):
        test_p3m_cpu_metallic = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.P3M,
            dict(prefactor=2., epsilon=0., mesh_off=[0.6, 0.7, 0.8], r_cut=1.5,
                 cao=2, mesh=[8, 10, 8], alpha=12., accuracy=0.01, tune=False,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12,
                 check_complex_residuals=False))
        test_p3m_cpu_non_metallic = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.P3M,
            dict(prefactor=2., epsilon=3., mesh_off=[0.6, 0.7, 0.8], r_cut=1.5,
                 cao=2, mesh=[8, 8, 8], alpha=12., accuracy=0.01, tune=False,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12))
        test_p3m_cpu_elc = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.ELC,
            dict(gap_size=2., maxPWerror=1e-3, const_pot=True, pot_diff=-3.,
                 delta_mid_top=0.5, delta_mid_bot=0.5, check_neutrality=False,
                 actor=espressomd.electrostatics.P3M(
                     prefactor=2., r_cut=1.5, cao=2, mesh=[8, 8, 8],
                     alpha=12., accuracy=0.01, tune=False)))

    if espressomd.has_features(["P3M", "CUDA"]) and espressomd.gpu_available():
        test_p3m_gpu_metallic = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.P3MGPU,
            dict(prefactor=2., epsilon=0., mesh_off=[0.6, 0.7, 0.8], r_cut=1.5,
                 cao=2, mesh=[8, 10, 8], alpha=12., accuracy=0.01, tune=False,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12,
                 check_complex_residuals=False))
        test_p3m_gpu_non_metallic = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.P3MGPU,
            dict(prefactor=2., epsilon=3., mesh_off=[0.6, 0.7, 0.8], r_cut=1.5,
                 cao=2, mesh=[8, 8, 8], alpha=12., accuracy=0.01, tune=False,
                 check_neutrality=True, charge_neutrality_tolerance=7e-12))
        test_p3m_gpu_elc = tests_common.generate_test_for_actor_class(
            system.electrostatics, espressomd.electrostatics.ELC,
            dict(gap_size=2., maxPWerror=1e-3, const_pot=True, pot_diff=-3.,
                 delta_mid_top=0.5, delta_mid_bot=0.5, check_neutrality=False,
                 actor=espressomd.electrostatics.P3MGPU(
                     prefactor=2., r_cut=1.5, cao=2, mesh=[8, 8, 8],
                     alpha=12., accuracy=0.01, tune=False)))

    def test_mmm1d_cpu(self):
        self.system.periodicity = [False, False, True]
        self.system.cell_system.set_n_square()
        valid_params = dict(
            prefactor=1., maxPWerror=1e-3, far_switch_radius=1.,
            check_neutrality=True, charge_neutrality_tolerance=7e-12,
            timings=5, verbose=False)
        tests_common.generate_test_for_actor_class(
            self.system.electrostatics, espressomd.electrostatics.MMM1D, valid_params)(self)

        for key in ["prefactor", "maxPWerror", "far_switch_radius", "timings"]:
            invalid_params = valid_params.copy()
            invalid_params[key] = -2
            with self.assertRaisesRegex(ValueError, f"Parameter '{key}' must be > 0"):
                espressomd.electrostatics.MMM1D(**invalid_params)

        # swapping two solvers should safely rollback to last valid solver
        self.assertEqual(abs(self.system.analysis.energy()["coulomb"]), 0.)
        mmm1d = espressomd.electrostatics.MMM1D(**valid_params)
        self.system.electrostatics.solver = mmm1d
        ref_energy = self.system.analysis.energy()["coulomb"]
        with self.assertRaisesRegex(RuntimeError, "CoulombP3M: requires periodicity"):
            p3m = espressomd.electrostatics.P3M(prefactor=1., accuracy=1e-2)
            self.system.electrostatics.solver = p3m
        self.assertEqual(self.system.electrostatics.solver, mmm1d)
        self.assertAlmostEqual(
            self.system.analysis.energy()["coulomb"], ref_energy, delta=1e-7)

    def test_charge_neutrality_check(self):
        self.system.part.add(pos=(0.0, 0.0, 0.0), q=1.)
        self.system.periodicity = [False, False, True]
        self.system.cell_system.set_n_square()
        actor = espressomd.electrostatics.MMM1D(prefactor=1.0, maxPWerror=1e-3)
        with self.assertRaisesRegex(RuntimeError, "The system is not charge neutral"):
            self.system.electrostatics.solver = actor
        self.assertIsNone(self.system.electrostatics.solver)
        self.assertFalse(actor.is_tuned)
        self.assertTrue(actor.check_neutrality)
        self.assertAlmostEqual(actor.charge_neutrality_tolerance, 2e-12)
        with self.assertRaisesRegex(ValueError, "Parameter 'charge_neutrality_tolerance' must be >= 0"):
            actor.charge_neutrality_tolerance = -1.
        actor.charge_neutrality_tolerance = None
        self.assertIsNone(actor.charge_neutrality_tolerance)
        self.assertFalse(actor.check_neutrality)
        actor.check_neutrality = True
        self.assertTrue(actor.check_neutrality)
        self.assertAlmostEqual(actor.charge_neutrality_tolerance, 2e-12)
        actor.check_neutrality = False
        self.assertFalse(actor.check_neutrality)
        self.assertIsNone(actor.charge_neutrality_tolerance)

    def test_mmm1d_cpu_tuning_exceptions(self):
        self.system.periodicity = [False, False, True]
        self.system.cell_system.set_n_square()
        actor = espressomd.electrostatics.MMM1D(
            prefactor=1., maxPWerror=1e-3, far_switch_radius=0.1)
        with self.assertRaisesRegex(RuntimeError, "MMM1D could not find a reasonable Bessel cutoff"):
            self.system.electrostatics.solver = actor
        self.assertIsNone(self.system.electrostatics.solver)
        self.assertFalse(actor.is_tuned)

    @utx.skipIfMissingFeatures(["P3M"])
    def test_elc_p3m_exceptions(self):
        P3M = espressomd.electrostatics.P3M
        ELC = espressomd.electrostatics.ELC
        P3MGPU = espressomd.electrostatics.P3MGPU
        # create valid solvers
        dh = espressomd.electrostatics.DH(prefactor=1.2, kappa=0.8, r_cut=2.0)
        p3m_params = dict(prefactor=1., epsilon=0.1, accuracy=1e-2,
                          mesh=[16, 16, 16], cao=7, r_cut=1.5, alpha=0.9)
        p3m = P3M(**p3m_params)
        elc = ELC(gap_size=2., maxPWerror=1., actor=p3m)

        # check runtime errors and input parameters
        with self.assertRaisesRegex(ValueError, "Parameter 'prefactor' must be > 0"):
            P3M(**{**p3m_params, 'prefactor': -2.})
        with self.assertRaisesRegex(ValueError, "Parameter 'timings' must be > 0"):
            P3M(**{**p3m_params, 'timings': -2})
        with self.assertRaisesRegex(ValueError, "Parameter 'mesh' has to be an integer or integer list of length 3"):
            P3M(**{**p3m_params, 'mesh': [8, 8]})
        if espressomd.has_features(["CUDA"]) and espressomd.gpu_available():
            with self.assertRaisesRegex(ValueError, "P3M GPU only implemented in single-precision mode"):
                P3MGPU(single_precision=False, **p3m_params)
        with self.assertRaisesRegex(ValueError, "Parameter 'actor' of type Coulomb::ElectrostaticLayerCorrection isn't supported by ELC"):
            ELC(gap_size=2., maxPWerror=1., actor=elc)
        with self.assertRaisesRegex(ValueError, "Parameter 'actor' of type Coulomb::DebyeHueckel isn't supported by ELC"):
            ELC(gap_size=2., maxPWerror=1., actor=dh)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'accuracy' is not a valid parameter"):
            ELC(gap_size=2., maxPWerror=1., actor=p3m, accuracy=1e-3)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'actor' is missing"):
            ELC(gap_size=2., maxPWerror=1.)
        with self.assertRaisesRegex(ValueError, "Parameter 'const_pot' must be True when 'pot_diff' is non-zero"):
            ELC(gap_size=2., maxPWerror=1., actor=p3m,
                const_pot=False, pot_diff=1.)
        with self.assertRaisesRegex(ValueError, "ELC with two parallel metallic boundaries requires the const_pot option"):
            ELC(gap_size=2., maxPWerror=1., actor=p3m, const_pot=False,
                delta_mid_top=-1., delta_mid_bot=-1.)
        # check contrasts are clamped if deviation from -1 is small
        epsilon = 1e-6
        with self.assertRaisesRegex(ValueError, r"Parameter 'delta_mid_top' must be >= -1 and <= \+1"):
            ELC(gap_size=2., maxPWerror=1., actor=p3m,
                delta_mid_top=-1. - epsilon, delta_mid_bot=-1.)
        with self.assertRaisesRegex(ValueError, r"Parameter 'delta_mid_bot' must be >= -1 and <= \+1"):
            ELC(gap_size=2., maxPWerror=1., actor=p3m,
                delta_mid_top=-1., delta_mid_bot=-1. - epsilon)
        epsilon = 1e-8
        elc = ELC(gap_size=2., maxPWerror=1., actor=p3m, const_pot=True,
                  delta_mid_top=-1. - epsilon, delta_mid_bot=1. + epsilon)
        self.assertAlmostEqual(elc.delta_mid_top, -1., delta=epsilon / 100.)
        self.assertAlmostEqual(elc.delta_mid_bot, +1., delta=epsilon / 100.)

        # run sanity checks
        elc = espressomd.electrostatics.ELC(
            gap_size=4., maxPWerror=1., actor=p3m)
        self.system.electrostatics.solver = elc
        periodicity_err_msg = r"requires periodicity \(True, True, True\)"
        with self.assertRaisesRegex(Exception, periodicity_err_msg):
            self.system.periodicity = [False, False, False]
        with self.assertRaisesRegex(Exception, periodicity_err_msg):
            self.system.integrator.run(0, recalc_forces=True)
        with self.assertRaisesRegex(Exception, periodicity_err_msg):
            self.system.analysis.energy()
        with self.assertRaisesRegex(Exception, periodicity_err_msg):
            self.system.analysis.pressure()
        self.system.periodicity = [True, True, True]
        n_nodes = self.system.cell_system.get_state()["n_nodes"]
        if n_nodes > 1:
            with self.assertRaisesRegex(Exception, "P3M: node grid must be sorted, largest first"):
                self.system.cell_system.node_grid = [1, n_nodes, 1]
            self.assertEqual(
                list(self.system.cell_system.node_grid),
                list(self.original_node_grid))
        with self.assertRaisesRegex(Exception, "ERROR: ELC gap size .+ larger than box length in z-direction"):
            self.system.change_volume_and_rescale_particles(2.5, "z")
        self.system.change_volume_and_rescale_particles(10., "z")
        self.system.electrostatics.solver = None
        with self.assertRaisesRegex(RuntimeError, "P3M real-space cutoff too large for ELC w/ dielectric contrast"):
            self.system.change_volume_and_rescale_particles(5., "z")
            elc = espressomd.electrostatics.ELC(
                actor=p3m,
                gap_size=1.,
                maxPWerror=1e-3,
                delta_mid_top=-1.,
                delta_mid_bot=-1.,
                const_pot=True,
                pot_diff=-3,
                check_neutrality=False,
            )
            self.system.electrostatics.solver = elc
        self.assertIsNone(self.system.electrostatics.solver)
        self.system.change_volume_and_rescale_particles(10., "z")
        self.system.periodicity = [True, True, False]
        with self.assertRaisesRegex(RuntimeError, periodicity_err_msg):
            elc = espressomd.electrostatics.ELC(
                gap_size=2., maxPWerror=1., actor=p3m)
            self.system.electrostatics.solver = elc
        self.assertIsNone(self.system.electrostatics.solver)
        self.system.periodicity = [True, True, True]


if __name__ == "__main__":
    ut.main()
