#
# Copyright (C) 2016-2022 The ESPResSo project
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
import espressomd
import numpy as np


@utx.skipIfMissingFeatures(["PARTICLE_ANISOTROPY", "MASS", "EXTERNAL_FORCES"])
class ThermostatsCommon(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    partcl_director = np.array(3 * [1 / np.sqrt(3.)])
    partcl_mass = 0.05

    gamma_parallel = 0.7
    gamma_ortho = 1.5
    kT = 0.123
    ext_force = 1.3 * np.array([0, 1, 0])

    n_steps_per_sample = 3
    n_samples = 10000
    seed = 42

    def tearDown(self) -> None:
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def setup_langevin(self, kT):
        self.system.thermostat.set_langevin(
            kT=kT,
            gamma=[
                self.gamma_ortho,
                self.gamma_ortho,
                self.gamma_parallel],
            seed=self.seed)
        self.system.integrator.set_vv()

    def setup_brownian(self, kT):
        self.system.thermostat.set_brownian(
            kT=kT,
            gamma=[
                self.gamma_ortho,
                self.gamma_ortho,
                self.gamma_parallel],
            seed=self.seed)
        self.system.integrator.set_brownian_dynamics()

    def setup_partcl(self, ext_force):
        return self.system.part.add(pos=self.system.box_l * np.random.random(3),
                                    director=self.partcl_director,
                                    mass=self.partcl_mass,
                                    rotation=3 * [False],
                                    ext_force=ext_force)

    def check_friction(self, partcl):
        vel = np.copy(partcl.v)
        force = np.copy(partcl.ext_force)
        force_parallel = np.dot(
            force, self.partcl_director) * self.partcl_director
        force_ortho = force - force_parallel
        vel_parallel = np.dot(vel, self.partcl_director) * self.partcl_director
        vel_ortho = vel - vel_parallel

        np.testing.assert_allclose(
            vel_parallel,
            force_parallel /
            self.gamma_parallel)
        np.testing.assert_allclose(vel_ortho, force_ortho / self.gamma_ortho)

    def test_friction_langevin(self):
        partcl = self.setup_partcl(self.ext_force)
        self.setup_langevin(0.)
        self.system.integrator.run(1000)
        self.check_friction(partcl)

    def test_friction_brownian(self):
        partcl = self.setup_partcl(self.ext_force)
        self.setup_brownian(0.)
        self.system.integrator.run(10)
        self.check_friction(partcl)

    def check_temperature(self, partcl):
        vels = []
        for _ in range(self.n_samples):
            self.system.integrator.run(self.n_steps_per_sample)
            vels.append(np.copy(partcl.v))
        vels = np.array(vels)

        vels_parallel = np.outer(
            np.dot(vels, self.partcl_director), self.partcl_director[None, :])
        vels_ortho = vels - vels_parallel
        vel_msq_parallel = np.mean(np.linalg.norm(vels_parallel, axis=1)**2)
        vel_msq_ortho = np.mean(np.linalg.norm(vels_ortho, axis=1)**2)

        # https://en.wikipedia.org/wiki/Thermal_velocity
        # parallel: 1d-diffusion
        kT_measured_parallel = self.partcl_mass * vel_msq_parallel
        # orthogonal: 2d-diffusion
        kT_measured_ortho = self.partcl_mass * vel_msq_ortho / 2

        for kT_measured in [kT_measured_parallel, kT_measured_ortho]:
            np.testing.assert_allclose(kT_measured, self.kT, rtol=3e-2)

    def test_temperature_brownian(self):
        partcl = self.setup_partcl(3 * [0.])
        self.setup_brownian(self.kT)
        self.check_temperature(partcl)

    def test_temperature_langevin(self):
        partcl = self.setup_partcl(3 * [0.])
        self.setup_langevin(self.kT)
        self.check_temperature(partcl)


if __name__ == "__main__":
    ut.main()
