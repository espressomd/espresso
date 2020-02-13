# Copyright (C) 2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import tests_common


@utx.skipIfMissingFeatures(["ROTATION", "PARTICLE_ANISOTROPY",
                            "ROTATIONAL_INERTIA", "DIPOLES"])
class RotDiffAniso(ut.TestCase):
    longMessage = True
    round_error_prec = 1E-14
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.skin = 5.0

    # The NVT thermostat parameters
    kT = 0.0
    gamma_global = np.zeros((3))

    # Particle properties
    J = [0.0, 0.0, 0.0]

    np.random.seed(4)

    def setUp(self):
        self.system.time = 0.0
        self.system.part.clear()
        if "BROWNIAN_DYNAMICS" in espressomd.features():
            self.system.thermostat.turn_off()
            # the default integrator is supposed implicitly
            self.system.integrator.set_nvt()

    def add_particles_setup(self, n):
        """
        Adding particles according to the
        previously set parameters.

        Parameters
        ----------
        n : :obj:`int`
            Number of particles.

        """

        for ind in range(n):
            part_pos = np.random.random(3) * self.box
            self.system.part.add(rotation=(1, 1, 1), id=ind,
                                 pos=part_pos)
            self.system.part[ind].rinertia = self.J
            if espressomd.has_features("ROTATION"):
                self.system.part[ind].omega_body = [0.0, 0.0, 0.0]

    def set_anisotropic_param(self):
        """
        Select parameters for anisotropic particles.

        """

        # NVT thermostat
        # Note: here & hereinafter specific variations in the random parameter
        # ranges are related to the test execution duration to achieve the
        # required statistical averages faster. The friction gamma_global should
        # be large enough in order to have the small enough D ~ kT / gamma and
        # to observe the details of the original rotational diffusion: the
        # Perrin1936 (see the reference below) tests are visible only when the
        # diffusive rotation is ~pi due to the exponential temporal dependencies
        # (see the equations referred in the check_rot_diffusion()).
        # Also, t0 ~ J / gamma should be small enough in order to remove the
        # small-time-scale diffusion effects which do not fit the Perrin1936's
        # tests which are based on the partial differential equation
        # (eq. (68), Perrin1934) leading only to the simple classical
        # Einstein-Smoluchowski equations of the diffusion in a contrast of the
        # eq. (10.2.26) [N. Pottier, doi:10.1007/s10955-010-0114-6 (2010)].
        self.gamma_global = 1E2 * np.random.uniform(0.35, 1.05, (3))

        # Particles' properties
        # As far as the problem characteristic time is t0 ~ J / gamma
        # and the Langevin equation finite-difference approximation is stable
        # only for time_step << t0, it is needed to set the moment of inertia
        # higher than some minimal value.
        # Also, it is expected to test the large enough J.
        # It should be not very large, otherwise the thermalization will require
        # too much of the CPU time: the in silico time should clock over the
        # t0.
        self.J = np.random.uniform(1.5, 16.5, (3))

    def set_isotropic_param(self):
        """
        Select parameters for isotropic particles.

        Parameters
        ----------

        """

        # NVT thermostat
        # see the comments in set_anisotropic_param()
        self.gamma_global[0] = 1E2 * np.random.uniform(0.35, 1.05)
        self.gamma_global[1] = self.gamma_global[0]
        self.gamma_global[2] = self.gamma_global[0]
        # Particles' properties
        # see the comments in set_anisotropic_param()
        self.J[0] = np.random.uniform(1.5, 16.5)
        self.J[1] = self.J[0]
        self.J[2] = self.J[0]

    def rot_diffusion_param_setup(self):
        """
        Setup the parameters for the rotational diffusion
        test check_rot_diffusion().

        """

        # Time
        # The time step should be less than t0 ~ mass / gamma
        self.system.time_step = 3E-3

        # Space
        self.box = 10.0
        self.system.box_l = 3 * [self.box]
        self.system.periodicity = [0, 0, 0]

        # NVT thermostat
        # Just some temperature range to cover by the test:
        self.kT = np.random.uniform(0.5, 1.5)

    def check_rot_diffusion(self, n):
        """
        The rotational diffusion tests based on the reference work
        [Perrin, F. (1936) Journal de Physique et Le Radium, 7(1), 1-11.
        https://doi.org/10.1051/jphysrad:01936007010100]
        with a theoretical background of
        [Perrin, F. (1934) Journal de Physique et Le Radium, 5(10), 497-511.
        https://doi.org/10.1051/jphysrad:01934005010049700]

        Parameters
        ----------
        n : :obj:`int`
            Number of particles.

        """
        # Global diffusivity tensor in the body frame:
        D = self.kT / self.gamma_global

        # Thermalizing...
        therm_steps = 100
        self.system.integrator.run(therm_steps)

        # Measuring...
        # Set the initial conditions according to the [Perrin1936], p.3.
        # The body angular velocity is rotated now, but there is only the
        # thermal velocity, hence, this does not impact the test and its
        # physical context.
        for ind in range(n):
            self.system.part[ind].quat = [1.0, 0.0, 0.0, 0.0]
        # Average direction cosines
        # Diagonal ones:
        dcosjj_validate = np.zeros((3))
        dcosjj_dev = np.zeros((3))
        # same to the power of 2
        dcosjj2_validate = np.zeros((3))
        dcosjj2_dev = np.zeros((3))
        # The non-diagonal elements for 2 different tests: negative ("nn") and
        # positive ("pp") ones.
        dcosijpp_validate = np.ones((3, 3))
        dcosijpp_dev = np.zeros((3, 3))
        dcosijnn_validate = np.ones((3, 3))
        dcosijnn_dev = np.zeros((3, 3))
        # The non-diagonal elements to the power of 2
        dcosij2_validate = np.ones((3, 3))
        dcosij2_dev = np.zeros((3, 3))

        self.system.time = 0.0
        int_steps = 20
        loops = 100
        for _ in range(loops):
            self.system.integrator.run(steps=int_steps)
            dcosjj = np.zeros((3))
            dcosjj2 = np.zeros((3))
            dcosijpp = np.zeros((3, 3))
            dcosijnn = np.zeros((3, 3))
            dcosij2 = np.zeros((3, 3))
            for ind in range(n):
                # Just a direction cosines functions averaging..
                dir_cos = tests_common.rotation_matrix_quat(self.system, ind)
                for j in range(3):
                    # the LHS of eq. (23) [Perrin1936].
                    dcosjj[j] += dir_cos[j, j]
                    # the LHS of eq. (32) [Perrin1936].
                    dcosjj2[j] += dir_cos[j, j]**2.0
                    for i in range(3):
                        if i != j:
                            # the LHS of eq. (24) [Perrin1936].
                            dcosijpp[i, j] += dir_cos[i, i] * dir_cos[j, j] + \
                                dir_cos[i, j] * dir_cos[j, i]
                            # the LHS of eq. (25) [Perrin1936].
                            dcosijnn[i, j] += dir_cos[i, i] * dir_cos[j, j] - \
                                dir_cos[i, j] * dir_cos[j, i]
                            # the LHS of eq. (33) [Perrin1936].
                            dcosij2[i, j] += dir_cos[i, j]**2.0
            dcosjj /= n
            dcosjj2 /= n
            dcosijpp /= n
            dcosijnn /= n
            dcosij2 /= n

            # Actual comparison.

            tolerance = 0.2
            # Too small values of the direction cosines are out of interest
            # compare to 0..1 range.
            min_value = 0.14

            # Eq. (23) [Perrin1936].
            dcosjj_validate[0] = np.exp(-(D[1] + D[2]) * self.system.time)
            dcosjj_validate[1] = np.exp(-(D[0] + D[2]) * self.system.time)
            dcosjj_validate[2] = np.exp(-(D[0] + D[1]) * self.system.time)
            dcosjj_dev = np.absolute(
                dcosjj - dcosjj_validate) / dcosjj_validate
            for j in range(3):
                if np.absolute(dcosjj_validate[j]) < min_value:
                    dcosjj_dev[j] = 0.0

            # Eq. (24) [Perrin1936].
            dcosijpp_validate[0, 1] = np.exp(
                -(4 * D[2] + D[1] + D[0]) * self.system.time)
            dcosijpp_validate[1, 0] = np.exp(
                -(4 * D[2] + D[1] + D[0]) * self.system.time)
            dcosijpp_validate[0, 2] = np.exp(
                -(4 * D[1] + D[2] + D[0]) * self.system.time)
            dcosijpp_validate[2, 0] = np.exp(
                -(4 * D[1] + D[2] + D[0]) * self.system.time)
            dcosijpp_validate[1, 2] = np.exp(
                -(4 * D[0] + D[2] + D[1]) * self.system.time)
            dcosijpp_validate[2, 1] = np.exp(
                -(4 * D[0] + D[2] + D[1]) * self.system.time)
            dcosijpp_dev = np.absolute(
                dcosijpp - dcosijpp_validate) / dcosijpp_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijpp_validate[i, j]) < min_value:
                        dcosijpp_dev[i, j] = 0.0

            # Eq. (25) [Perrin1936].
            dcosijnn_validate[0, 1] = np.exp(-(D[1] + D[0]) * self.system.time)
            dcosijnn_validate[1, 0] = np.exp(-(D[1] + D[0]) * self.system.time)
            dcosijnn_validate[0, 2] = np.exp(-(D[2] + D[0]) * self.system.time)
            dcosijnn_validate[2, 0] = np.exp(-(D[2] + D[0]) * self.system.time)
            dcosijnn_validate[1, 2] = np.exp(-(D[2] + D[1]) * self.system.time)
            dcosijnn_validate[2, 1] = np.exp(-(D[2] + D[1]) * self.system.time)
            dcosijnn_dev = np.absolute(
                dcosijnn - dcosijnn_validate) / dcosijnn_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosijnn_validate[i, j]) < min_value:
                        dcosijnn_dev[i, j] = 0.0

            # Eq. (30) [Perrin1936].
            D0 = sum(D[:]) / 3.0
            D1D1 = 0.0
            for j in range(3):
                for i in range(3):
                    if i != j:
                        D1D1 += D[i] * D[j]
            D1D1 /= 6.0
            # Technical workaround of a digital arithmetic issue for isotropic
            # particle
            if np.absolute((D0**2 - D1D1) / (D0**2 + D1D1)
                           ) < self.round_error_prec:
                D1D1 *= (1.0 - 2.0 * self.round_error_prec)
            # Eq. (32) [Perrin1936].
            dcosjj2_validate = 1. / 3. + (1. / 3.) * (1. + (D - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                + (1. / 3.) * (1. - (D - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosjj2_dev = np.absolute(
                dcosjj2 - dcosjj2_validate) / dcosjj2_validate
            for j in range(3):
                if np.absolute(dcosjj2_validate[j]) < min_value:
                    dcosjj2_dev[j] = 0.0

            # Eq. (33) [Perrin1936].
            dcosij2_validate[0, 1] = 1. / 3. - (1. / 6.) * (1. - (D[2] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[2] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[1, 0] = dcosij2_validate[0, 1]

            dcosij2_validate[0, 2] = 1. / 3. - (1. / 6.) * (1. - (D[1] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[1] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[2, 0] = dcosij2_validate[0, 2]

            dcosij2_validate[1, 2] = 1. / 3. - (1. / 6.) * (1. - (D[0] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 - np.sqrt(D0**2 - D1D1)) * self.system.time) \
                - (1. / 6.) * (1. + (D[0] - D0) / (2. * np.sqrt(D0**2 - D1D1))) \
                * np.exp(-6. * (D0 + np.sqrt(D0**2 - D1D1)) * self.system.time)
            dcosij2_validate[2, 1] = dcosij2_validate[1, 2]
            dcosij2_dev = np.absolute(
                dcosij2 - dcosij2_validate) / dcosij2_validate
            for i in range(3):
                for j in range(3):
                    if np.absolute(dcosij2_validate[i, j]) < min_value:
                        dcosij2_dev[i, j] = 0.0

            for j in range(3):
                self.assertLessEqual(
                    abs(dcosjj_dev[j]), tolerance,
                    msg='Relative deviation dcosjj_dev[{0}] in a rotational '
                        'diffusion is too large: {1}'.format(j, dcosjj_dev[j]))
                self.assertLessEqual(
                    abs(dcosjj2_dev[j]), tolerance,
                    msg='Relative deviation dcosjj2_dev[{0}] in a rotational '
                        'diffusion is too large: {1}'.format(j, dcosjj2_dev[j]))
                for i in range(3):
                    if i != j:
                        self.assertLessEqual(
                            abs(dcosijpp_dev[i, j]), tolerance,
                            msg='Relative deviation dcosijpp_dev[{0},{1}] in a '
                                'rotational diffusion is too large: {2}'
                                .format(i, j, dcosijpp_dev[i, j]))
                        self.assertLessEqual(
                            abs(dcosijnn_dev[i, j]), tolerance,
                            msg='Relative deviation dcosijnn_dev[{0},{1}] in a '
                                'rotational diffusion is too large: {2}'
                                .format(i, j, dcosijnn_dev[i, j]))
                        self.assertLessEqual(
                            abs(dcosij2_dev[i, j]), tolerance,
                            msg='Relative deviation dcosij2_dev[{0},{1}] in a '
                                'rotational diffusion is too large: {2}'
                                .format(i, j, dcosij2_dev[i, j]))

    # Langevin Dynamics / Anisotropic
    def test_case_00(self):
        n = 800
        self.rot_diffusion_param_setup()
        self.set_anisotropic_param()
        self.add_particles_setup(n)
        self.system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        # Actual integration and validation run
        self.check_rot_diffusion(n)

    # Langevin Dynamics / Isotropic
    def test_case_01(self):
        n = 800
        self.rot_diffusion_param_setup()
        self.set_isotropic_param()
        self.add_particles_setup(n)
        self.system.thermostat.set_langevin(
            kT=self.kT, gamma=self.gamma_global, seed=42)
        # Actual integration and validation run
        self.check_rot_diffusion(n)

    if "BROWNIAN_DYNAMICS" in espressomd.features():
        # Brownian Dynamics / Isotropic
        def test_case_10(self):
            n = 800
            self.system.thermostat.turn_off()
            self.rot_diffusion_param_setup()
            self.set_isotropic_param()
            self.add_particles_setup(n)
            self.system.thermostat.set_brownian(
                kT=self.kT, gamma=self.gamma_global, seed=42)
            self.system.integrator.set_brownian_dynamics()
            # Actual integration and validation run
            self.check_rot_diffusion(n)


if __name__ == '__main__':
    ut.main()
