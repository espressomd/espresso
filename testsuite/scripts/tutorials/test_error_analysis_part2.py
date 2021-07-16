# Copyright (C) 2021 The ESPResSo project
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
import importlib_wrapper
import numpy as np
import scipy.signal

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    filepath="@TUTORIALS_DIR@/error_analysis/error_analysis_part2.py",
    gpu=False,
    random_seeds=True)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def ar_1_process(self, n, c, phi, eps):
        y0 = np.random.normal(loc=c / (1 - phi),
                              scale=np.sqrt(eps**2 / (1 - phi**2)))
        y = c + np.random.normal(loc=0.0, scale=eps, size=n - 1)
        y = np.insert(y, 0, y0)
        # get an AR(1) process from an ARMA(p,q) process with p=1 and q=0
        y = scipy.signal.lfilter([1.], [1., -phi], y)
        return y

    def test_ar1_implementation(self):
        with self.assertRaises(ValueError):
            tutorial.ar_1_process(10, 1.0, 1.1, 3.0)
        with self.assertRaises(ValueError):
            tutorial.ar_1_process(10, 1.0, -1.1, 3.0)

        for seed in range(5):
            for eps in [0.5, 1., 2.]:
                for phi in [0.1, 0.8, 0.999, -0.3]:
                    c = eps / 2.
                    np.random.seed(seed)
                    seq = tutorial.ar_1_process(10, c, phi, eps)
                    np.random.seed(seed)
                    ref = self.ar_1_process(10, c, phi, eps)
                    np.testing.assert_allclose(seq, ref, atol=1e-12, rtol=0)

    def test(self):
        self.assertLess(abs(tutorial.PHI_1), 1.0)
        self.assertLess(abs(tutorial.PHI_2), 1.0)

        # The analytic expressions for the AR(1) process are taken from
        # https://en.wikipedia.org/wiki/Autoregressive_model#Example:_An_AR(1)_process
        # (accessed June 2021)
        SIGMA_1 = np.sqrt(tutorial.EPS_1 ** 2 / (1 - tutorial.PHI_1 ** 2))
        TAU_EXP_1 = -1 / np.log(tutorial.PHI_1)

        np.testing.assert_allclose(tutorial.an_acf_1,
                                   SIGMA_1**2 * np.exp(-np.linspace(0,
                                                                    tutorial.N_MAX - 1,
                                                                    tutorial.N_MAX) / TAU_EXP_1))
        # The autocorrelation is exponential, thus tau_exp = tau_int, and
        # therefore
        N_EFF_1 = tutorial.N_SAMPLES / (2 * TAU_EXP_1)
        SEM_1 = np.sqrt(SIGMA_1 ** 2 / N_EFF_1)

        self.assertAlmostEqual(tutorial.sem, SEM_1, delta=0.1 * SEM_1)
        self.assertAlmostEqual(tutorial.N_eff, N_EFF_1, delta=0.1 * N_EFF_1)
        # for some reason, the integrated autocorrelation time is always higher
        # than the exponential one, in the tutorial
        self.assertAlmostEqual(
            tutorial.tau_int,
            TAU_EXP_1,
            delta=0.1 * TAU_EXP_1)

        SIGMA_2 = np.sqrt(tutorial.EPS_2 ** 2 / (1 - tutorial.PHI_2 ** 2))
        TAU_EXP_2 = -1 / np.log(tutorial.PHI_2)
        SEM_2 = np.sqrt(2 * SIGMA_2 ** 2 * TAU_EXP_2 / tutorial.N_SAMPLES)
        # the point of the following value in the tutorial is that it is very
        # inaccurate, thus the high tolerance
        self.assertAlmostEqual(tutorial.sem_2, SEM_2, delta=0.2 * SEM_2)


if __name__ == "__main__":
    ut.main()
