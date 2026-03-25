#
# Copyright (C) 2020-2026 The ESPResSo project
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
import espressomd.interactions
import espressomd.lees_edwards
import espressomd.shapes
import espressomd.propagation
import os
import numpy as np
import unittest as ut
import unittest_decorators as utx


class Test(ut.TestCase):

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.01
    msg = r'while calling method integrate\(\): ERROR: '

    def setUp(self):
        self.system.box_l = [1., 1., 1.]
        self.system.part.add(pos=(0, 0, 0))
        self.system.integrator.set_vv()
        self.system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()
        self.system.constraints.clear()
        self.system.lees_edwards.protocol = None

    def test_00_common_interface(self):
        Propagation = espressomd.propagation.Propagation
        self.system.integrator.set_vv()
        with self.assertRaisesRegex(ValueError, 'time_step must be > 0.'):
            self.system.time_step = -2.
        with self.assertRaisesRegex(ValueError, 'time_step must be > 0.'):
            self.system.time_step = 0.
        with self.assertRaisesRegex(ValueError, "Parameter 'steps' must be positive"):
            self.system.integrator.run(steps=-1)
        with self.assertRaisesRegex(RuntimeError, 'cannot automatically determine skin, please set it manually'):
            self.system.integrator.run()
        self.system.cell_system.skin = 0.4
        with self.assertRaisesRegex(ValueError, 'cannot reuse old forces and recalculate forces'):
            self.system.integrator.run(recalc_forces=True, reuse_forces=True)
        self.assertIsNone(self.system.integrator.integrator.call_method("unk"))
        self.assertIsNone(self.system.thermostat.call_method("unk"))
        self.assertIsNone(self.system.thermostat.kT)
        if espressomd.has_features("WCA"):
            wca = self.system.non_bonded_inter[0, 0].wca
            wca.set_params(epsilon=1., sigma=0.01)
            wall = espressomd.shapes.Wall(normal=[0., 0., 1.], dist=100.)
            self.system.constraints.add(
                shape=wall, particle_type=0, penetrable=False)
            with self.assertRaisesRegex(Exception, self.msg + 'Constraint violated by particle 0 dist -100'):
                self.system.integrator.run(0)
            with self.assertRaisesRegex(Exception, 'Constraint violated by particle 0'):
                self.system.analysis.energy()
            wca.set_params(epsilon=0., sigma=0.)
        if espressomd.has_features(["ROTATION"]):
            p = self.system.part.by_id(0)
            with self.assertRaisesRegex(Exception, "Rotating particles must have a rotation propagation mode enabled"):
                p.propagation = Propagation.TRANS_LANGEVIN
                p.rotation = [False, False, True]
                self.system.integrator.run(0, recalc_forces=True)
            p.propagation = Propagation.SYSTEM_DEFAULT
        if espressomd.has_features(["THERMAL_STONER_WOHLFARTH"]):
            p0 = self.system.part.by_id(0)
            p1 = self.system.part.add(pos=p0.pos)
            p1.vs_auto_relate_to(p0)
            p1.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT
            magnetodynamics = p1.magnetodynamics
            with self.assertRaisesRegex(Exception, "The thermal Stoner-Wohlfarth model requires the Langevin thermostat"):
                magnetodynamics["is_enabled"] = True
                p1.magnetodynamics = magnetodynamics
                self.system.integrator.run(0, recalc_forces=True)
            magnetodynamics["is_enabled"] = False
            p1.magnetodynamics = magnetodynamics
            self.system.integrator.run(0, recalc_forces=True)

    def test_01_statefulness(self):
        # setting a thermostat with invalid values should be a no-op
        self.assertIsNone(self.system.thermostat.kT)
        self.assertIsNone(self.system.thermostat.langevin.seed)
        self.assertIsNone(self.system.thermostat.langevin.gamma)
        with self.assertRaisesRegex(ValueError, "Parameter 'seed' must be a positive integer"):
            self.system.thermostat.set_langevin(kT=1., gamma=1., seed=-1)
        self.assertIsNone(self.system.thermostat.kT)
        self.assertIsNone(self.system.thermostat.langevin.seed)
        self.assertIsNone(self.system.thermostat.langevin.gamma)
        with self.assertRaisesRegex(ValueError, "Parameter 'kT' cannot be negative"):
            self.system.thermostat.set_langevin(kT=-1., gamma=1., seed=42)
        self.assertIsNone(self.system.thermostat.kT)
        self.assertIsNone(self.system.thermostat.langevin.seed)
        self.assertIsNone(self.system.thermostat.langevin.gamma)
        with self.assertRaisesRegex(ValueError, "Parameter 'gamma' cannot be negative"):
            self.system.thermostat.set_langevin(kT=1., gamma=-1., seed=42)
        self.assertIsNone(self.system.thermostat.kT)
        self.assertIsNone(self.system.thermostat.langevin.seed)
        self.assertIsNone(self.system.thermostat.langevin.gamma)

    def test_vv_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_vv()
        with self.assertRaisesRegex(Exception, self.msg + 'The VV integrator is incompatible with the currently active combination of thermostats'):
            self.system.integrator.run(0)

    def test_brownian_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.integrator.set_brownian_dynamics()
        self.assertIsNone(self.system.thermostat.kT)
        with self.assertRaisesRegex(Exception, self.msg + 'The BD integrator requires the BD thermostat'):
            self.system.integrator.run(0)
        with self.assertRaisesRegex(RuntimeError, "Parameter 'brownian' is read-only."):
            self.system.thermostat.brownian = 1
        self.assertIsNone(self.system.thermostat.kT)

    def test_langevin_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.integrator.set_vv()
        self.system.thermostat.set_langevin(kT=2., gamma=3., seed=42)

        def check_original_params(thermo_off):
            langevin = self.system.thermostat.langevin
            np.testing.assert_allclose(np.copy(langevin.gamma), 3.)
            np.testing.assert_allclose(np.copy(langevin.seed), 42)
            if thermo_off:
                self.assertIsNone(self.system.thermostat.kT)
            else:
                np.testing.assert_allclose(self.system.thermostat.kT, 2.)

        # updating a thermostat with invalid parameters should raise an
        # exception and roll back to the last valid state of the thermostat
        for thermo_off in [False, True]:
            if thermo_off:
                self.system.thermostat.turn_off()
            with self.assertRaisesRegex(ValueError, "Parameter 'seed' must be a positive integer"):
                self.system.thermostat.set_langevin(kT=1., gamma=1., seed=-1)
            check_original_params(thermo_off)
            with self.assertRaisesRegex(ValueError, "Parameter 'kT' cannot be negative"):
                self.system.thermostat.set_langevin(kT=-1., gamma=1., seed=42)
            check_original_params(thermo_off)
            with self.assertRaisesRegex(ValueError, "Parameter 'gamma' cannot .* negative"):
                self.system.thermostat.set_langevin(kT=1., gamma=-1., seed=42)
            check_original_params(thermo_off)

        with self.assertRaisesRegex(RuntimeError, "Parameter 'langevin' is read-only."):
            self.system.thermostat.langevin = 1
        with self.assertRaisesRegex(RuntimeError, "Parameter 'gamma' is read-only."):
            self.system.thermostat.langevin.gamma = 1
        with self.assertRaisesRegex(RuntimeError, "Parameter 'seed' is read-only."):
            self.system.thermostat.langevin.seed = 1
        with self.assertRaisesRegex(RuntimeError, "Parameter 'philox_counter' is read-only."):
            self.system.thermostat.langevin.philox_counter = 1

    @utx.skipIfMissingFeatures("NPT")
    def test_npt_integrator(self):
        self.system.cell_system.skin = 0.4
        with self.assertRaisesRegex(Exception, "Parameter 'barostat' must be 'Andersen' or 'MTK'."):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1.0, piston=1.0, barostat='XXX')
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)
        with self.assertRaisesRegex(Exception, self.msg + 'The NpT integrator requires the NpT thermostat'):
            self.system.integrator.run(0)
        self.system.thermostat.turn_off()
        self.system.thermostat.set_npt(kT=1.0, gamma0=2., gammav=0.04, seed=42)
        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=1., initial_pos_offset=0., time_0=0.)
        with self.assertRaisesRegex(Exception, self.msg + 'The NpT integrator cannot use Lees-Edwards'):
            self.system.lees_edwards.set_boundary_conditions(
                shear_direction="x", shear_plane_normal="y", protocol=protocol)
            self.system.integrator.run(0)
        with self.assertRaisesRegex(Exception, self.msg + 'The NpT integrator cannot use Lees-Edwards'):
            self.system.lees_edwards.protocol = espressomd.lees_edwards.LinearShear(
                shear_velocity=0., initial_pos_offset=0., time_0=0.)
            self.system.integrator.run(0)
        self.system.lees_edwards.protocol = None
        self.system.integrator.run(0)

    @ut.skipIf(espressomd.code_info._CodeInfo().call_method("has_fast_math"),
               "cannot run with fast-math optimizations")
    @ut.skipIf(os.environ.get("UBSAN_OPTIONS"),
               "cannot run with UBSAN instrumentation")
    @ut.skipIf(espressomd.has_features("FPE"),
               "cannot run with FPE instrumentation")
    @utx.skipIfMissingFeatures(["NPT", "WCA"])
    def test_npt_integrator_negative_volume(self):
        """Test for NpT with bad parameters."""

        import tests_common
        data = np.genfromtxt(tests_common.data_path("npt_lj_system.data"))
        ref_box_l = np.max(data[:, 0:3])

        system = self.system
        system.part.clear()
        system.cell_system.skin = 0.

        for barostat in ["Andersen", "MTK"]:
            system.box_l = 3 * [ref_box_l]
            system.time_step = 0.01
            if barostat == "Andersen":
                piston = 1e-4
            else:
                piston = 4.0
            direction = [True] * 3
            ext_pressure = 100.0  # Too large external pressure
            system.part.add(pos=data[:, 0:3], v=data[:, 3:6])
            system.integrator.set_vv()
            system.thermostat.set_npt(kT=1.0, gamma0=0.1, gammav=1e-3, seed=42)
            system.integrator.set_isotropic_npt(ext_pressure=ext_pressure,
                                                piston=piston,
                                                direction=direction,
                                                barostat=barostat)

            if barostat == "Andersen":
                with self.assertRaisesRegex(Exception, "the volume to become negative"):
                    system.integrator.run(10)
                system.part.clear()
            if barostat == "MTK":
                # Volume cannot be negative with NPT based on MTK equation
                self.assertGreater(float(np.prod(system.box_l)), 0.)

        system.part.clear()
        system.box_l = [1., 1., 1.]

    @utx.skipIfMissingFeatures("STOKESIAN_DYNAMICS")
    def test_stokesian_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.periodicity = 3 * [False]
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_stokesian_dynamics(
            viscosity=1.0, radii={0: 1.0})
        with self.assertRaisesRegex(Exception, self.msg + 'The SD integrator requires the SD thermostat'):
            self.system.integrator.run(0)

    def test_steepest_descent_integrator(self):
        self.system.cell_system.skin = 0.4
        params = {"f_max": 0., "gamma": 0.1, "max_displacement": 5.}
        for key in params:
            invalid_params = params.copy()
            del invalid_params[key]
            with self.assertRaisesRegex(RuntimeError, f"Parameter '{key}' is missing"):
                self.system.integrator.set_steepest_descent(
                    **invalid_params)
        self.system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        with self.assertRaisesRegex(Exception, self.msg + 'The steepest descent integrator is incompatible with thermostats'):
            self.system.integrator.run(0)

    def test_temperature_change(self):
        # temperature change only allowed when no other thermostat is active
        self.system.thermostat.set_langevin(kT=1., gamma=1., seed=42)
        self.system.thermostat.set_langevin(kT=2., gamma=1., seed=42)
        self.system.thermostat.set_brownian(kT=2., gamma=1., seed=42)
        with self.assertRaisesRegex(RuntimeError, "Cannot set parameter 'kT' to 1.0*: there are currently active thermostats with kT=2.0*"):
            self.system.thermostat.set_brownian(kT=1., gamma=1., seed=42)
        with self.assertRaisesRegex(RuntimeError, f"Parameter 'kT' is read-only"):
            self.system.thermostat.kT = 2.


if __name__ == "__main__":
    ut.main()
