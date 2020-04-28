#
# Copyright (C) 2013-2019 The ESPResSo project
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
import numpy as np
import sys
import unittest as ut
import unittest_decorators as utx
import tests_common
import espressomd


if not espressomd.has_features("ELECTROSTATICS"):
    sys.exit()


class ElectrostaticInteractionsTests:
    # Handle to espresso system
    system = espressomd.System(box_l=[10.0] * 3)
    system.periodicity = [0, 0, 1]
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.cell_system.set_n_square()
    system.thermostat.set_langevin(kT=0, gamma=1, seed=8)

    pid_target, pos_x_target, pos_y_target, pos_z_target, q_target, f_x_target, f_y_target, f_z_target = np.loadtxt(
        tests_common.abspath("data/mmm1d_data.txt"), unpack=True)
    vec_f_target = np.stack((f_x_target, f_y_target, f_z_target), axis=-1)
    energy_target = -7.156365298205383
    num_particles = pid_target.shape[0]

    allowed_error = 1e-4

    def setUp(self):
        for i in range(self.num_particles):
            self.system.part.add(
                pos=[self.pos_x_target[i],
                     self.pos_y_target[i],
                     self.pos_z_target[i]],
                q=self.q_target[i])
        self.mmm1d = self.MMM1D(prefactor=1.0, maxPWerror=1e-20)
        self.system.actors.add(self.mmm1d)
        self.system.integrator.run(steps=0)

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def test_forces(self):
        measured_f = self.system.part[:].f
        for i in range(self.num_particles):
            for comp in range(3):
                self.assertAlmostEqual(
                    measured_f[i, comp], self.vec_f_target[i, comp],
                    delta=self.allowed_error,
                    msg="Measured force deviates too much "
                        "for particle {} in component {}".format(i, comp))

    def test_energy(self):
        measured_el_energy = self.system.analysis.energy()["total"] \
            - self.system.analysis.energy()["kinetic"]
        self.assertAlmostEqual(
            measured_el_energy, self.energy_target, delta=self.allowed_error,
            msg="Measured energy deviates too much from stored result")

    def test_with_analytical_result(self, prefactor=1.0, accuracy=1e-4):
        self.system.part.clear()
        self.system.part.add(pos=[0, 0, 0], q=1)
        self.system.part.add(pos=[0, 0, 1], q=1)

        self.system.integrator.run(steps=0)
        f_measured = self.system.part[0].f
        energy_measured = self.system.analysis.energy()["total"]
        target_energy_config = 1.00242505606 * prefactor
        target_force_z_config = -0.99510759 * prefactor

        self.assertAlmostEqual(
            f_measured[0], 0, delta=self.allowed_error,
            msg="Measured force in x deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[1], 0, delta=self.allowed_error,
            msg="Measured force in y deviates too much from analytical result")
        self.assertAlmostEqual(
            f_measured[2], target_force_z_config, delta=accuracy,
            msg="Measured force in z deviates too much from analytical result")
        self.assertAlmostEqual(
            energy_measured, target_energy_config, delta=self.allowed_error,
            msg="Measured energy deviates too much from analytical result")

    def test_bjerrum_length_change(self):
        self.system.part.clear()
        self.system.actors.clear()
        prefactor = 2
        mmm1d = self.MMM1D(prefactor=prefactor, maxPWerror=1e-20)
        self.system.actors.add(mmm1d)
        self.test_with_analytical_result(prefactor=prefactor, accuracy=0.0017)


@utx.skipIfMissingFeatures(["ELECTROSTATICS", "MMM1D_GPU"])
class MMM1D_GPU_Test(ElectrostaticInteractionsTests, ut.TestCase):
    from espressomd.electrostatics import MMM1D


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class MMM1D_Test(ElectrostaticInteractionsTests, ut.TestCase):
    from espressomd.electrostatics import MMM1D


if __name__ == "__main__":
    ut.main()
