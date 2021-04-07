#
# Copyright (C) 2019-2020 The ESPResSo project
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
import espressomd.observables
import espressomd.accumulators
import espressomd.constraints
import numpy as np
import unittest as ut
import unittest_decorators as utx
from tests_common import abspath

s = espressomd.System(box_l=[1.0, 1.0, 1.0])


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsTest(ut.TestCase):
    system = s

    # Digitized reference data of Figure 5b from
    # Durlofsky et al., J. Fluid Mech. 180, 21 (1987)
    # https://doi.org/10.1017/S002211208700171X
    data = np.loadtxt(abspath('data/dancing.txt'))

    def setUp(self):
        self.system.box_l = [10] * 3
        self.system.periodicity = [0, 0, 0]
        self.system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()

    def falling_spheres(self, time_step, l_factor, t_factor, sd_method='fts'):
        self.system.time_step = time_step
        self.system.part.add(pos=[-5 * l_factor, 0, 0], rotation=[1, 1, 1])
        self.system.part.add(pos=[0 * l_factor, 0, 0], rotation=[1, 1, 1])
        self.system.part.add(pos=[7 * l_factor, 0, 0], rotation=[1, 1, 1])

        self.system.integrator.set_stokesian_dynamics(
            viscosity=1.0 / (t_factor * l_factor),
            radii={0: 1.0 * l_factor}, approximation_method=sd_method)

        gravity = espressomd.constraints.Gravity(
            g=[0, -1.0 * l_factor / (t_factor**2), 0])
        self.system.constraints.add(gravity)
        self.system.time_step = 1.0 * t_factor

        obs = espressomd.observables.ParticlePositions(ids=(0, 1, 2))
        acc = espressomd.accumulators.TimeSeries(obs=obs, delta_N=1)
        self.system.auto_update_accumulators.add(acc)
        acc.update()

        if sd_method == 'fts':
            y_min = -555
            intsteps = 8000
        else:
            y_min = -200
            intsteps = 3000
        intsteps = int(intsteps * t_factor / self.system.time_step)

        self.system.integrator.run(intsteps)

        simul = acc.time_series()[:, :, 0:2]
        paper = self.data.reshape([-1, 3, 2])

        for pid in range(3):
            dist = []
            # the simulated trajectory is oversampled by a ratio of
            # (90/t_factor):1 compared to the published trajectory
            for desired in paper[:, pid] * l_factor:
                if desired[1] < y_min * l_factor:
                    break
                # find the closest point in the simulated trajectory
                idx = np.abs(simul[:, pid, 1] - desired[1]).argmin()
                actual = simul[idx, pid]
                dist.append(np.linalg.norm(actual - desired))
            self.assertLess(idx, intsteps, msg='Insufficient sampling')
            np.testing.assert_allclose(dist, 0, rtol=0, atol=0.5 * l_factor)

    def test_default(self):
        self.falling_spheres(1.0, 1.0, 1.0)

    def test_rescaled(self):
        self.falling_spheres(1.0, 4.5, 2.5)

    def test_different_time_step(self):
        self.falling_spheres(0.7, 1.0, 1.0)

    def test_default_ft(self):
        self.falling_spheres(1.0, 1.0, 1.0, 'ft')


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDiffusionTest(ut.TestCase):
    system = s

    kT = 1e-4
    R = 1.5
    eta = 2.4

    def setUp(self):
        self.system.box_l = [10] * 3
        self.system.periodicity = [0, 0, 0]

        self.system.time_step = 1.0
        self.system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()
        self.system.thermostat.set_stokesian(kT=0)

    def test(self):
        p = self.system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
        self.system.integrator.set_stokesian_dynamics(
            viscosity=self.eta, radii={0: self.R})
        self.system.thermostat.set_stokesian(kT=self.kT, seed=42)

        intsteps = int(100000 / self.system.time_step)
        pos = np.empty([intsteps + 1, 3])
        orientation = np.empty((intsteps + 1, 3))

        pos[0, :] = p.pos
        orientation[0, :] = p.director
        for i in range(intsteps):
            self.system.integrator.run(1)
            pos[i + 1, :] = p.pos
            orientation[i + 1, :] = p.director

        # translational diffusion coefficient
        D_expected = self.kT / (6 * np.pi * self.eta * self.R)

        # NOTE on steps_per_slice:
        # The shorter these trajectories are, the more trajectories we get
        # and the more accurate the result is.
        # However since we want to test that diffusion works for
        # "long" trajectories we won't go too low.
        n_steps_per_slice = 200
        n_slices = int(intsteps / n_steps_per_slice)

        squared_displacement_per_slice = np.empty(n_slices)
        for i in range(n_slices):
            squared_displacement_per_slice[i] = np.linalg.norm(
                pos[i * n_steps_per_slice] - pos[(i + 1) * n_steps_per_slice])**2
        D_measured = np.mean(squared_displacement_per_slice) / \
            (6 * n_steps_per_slice * self.system.time_step)
        self.assertAlmostEqual(D_expected, D_measured, delta=D_expected * 0.1)

        # rotational diffusion coefficient
        Dr_expected = self.kT / (8 * np.pi * self.eta * self.R**3)
        n_steps_per_slice = 200
        n_slices = int(intsteps / n_steps_per_slice)

        # This is the displacement equivalent to measure rotational diffusion
        costheta_per_slice = np.empty(n_slices)

        for i in range(n_slices):
            costheta_per_slice[i] = np.dot(orientation[i * n_steps_per_slice, :],
                                           orientation[(i + 1) * n_steps_per_slice, :])
        Dr_measured = -np.log(np.mean(costheta_per_slice)) / \
            (2 * n_steps_per_slice * self.system.time_step)
        self.assertAlmostEqual(
            Dr_expected,
            Dr_measured,
            delta=Dr_expected * 0.1)


if __name__ == '__main__':
    ut.main()
