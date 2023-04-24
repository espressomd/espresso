#
# Copyright (C) 2020-2022 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx


class Test(ut.TestCase):

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.01
    msg = r'while calling method integrate\(\): ERROR: '

    def setUp(self):
        self.system.part.add(pos=(0, 0, 0))
        self.system.integrator.set_vv()
        self.system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.part.clear()
        self.system.constraints.clear()
        self.system.lees_edwards.protocol = None

    def test_00_common_interface(self):
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

    def test_vv_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=42)
        self.system.integrator.set_vv()
        with self.assertRaisesRegex(Exception, self.msg + 'The VV integrator is incompatible with the currently active combination of thermostats'):
            self.system.integrator.run(0)

    def test_brownian_integrator(self):
        self.system.cell_system.skin = 0.4
        self.system.integrator.set_brownian_dynamics()
        with self.assertRaisesRegex(Exception, self.msg + 'The BD integrator requires the BD thermostat'):
            self.system.integrator.run(0)

    @utx.skipIfMissingFeatures("NPT")
    def test_npt_integrator(self):
        self.system.cell_system.skin = 0.4
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
            self.system.lees_edwards.protocol = espressomd.lees_edwards.Off()
            self.system.integrator.run(0)
        self.system.lees_edwards.protocol = None
        self.system.integrator.run(0)

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


if __name__ == "__main__":
    ut.main()
