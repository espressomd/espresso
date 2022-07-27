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

import unittest as ut
import numpy as np
import espressomd.analyze


class AnalyzeAcf(ut.TestCase):
    np.random.seed(42)

    def test_acf(self):
        # the ACF of white noise is a Dirac function
        N = 2**14
        white_noise = np.random.uniform(-1, 1, N)
        acf_ref = np.zeros(N)
        acf_ref[0] = 1
        acf = espressomd.analyze.autocorrelation(white_noise)
        acf_normalized = acf / np.var(white_noise) * (N - np.arange(N)) / N
        np.testing.assert_allclose(acf_normalized, acf_ref, atol=0.05, rtol=0.)
        # the ACF of an auto-regressive model AR(p) of order p is a sum of
        # decaying exponentials, ACF(AR(1)) is a single exponential
        phi = -0.9
        ar1 = np.copy(white_noise)
        ar1[0] = 200
        for i in range(1, N):
            ar1[i] += phi * ar1[i - 1]
        acf_ref = phi**np.arange(N)[:50]
        acf = espressomd.analyze.autocorrelation(ar1)
        acf_normalized = acf[:50] / acf[0]
        np.testing.assert_allclose(acf_normalized, acf_ref, atol=0.05, rtol=0.)
        # the ACF of 2-dimensional data is the sum of 1-dimensional ACFs
        data_2d = np.array([np.random.uniform(-1, 1, N),
                            np.random.uniform(-1, 1, N),
                            np.random.uniform(-1, 1, N)]).T
        acf = espressomd.analyze.autocorrelation(data_2d)
        acf_ref = (espressomd.analyze.autocorrelation(data_2d[:, 0]) +
                   espressomd.analyze.autocorrelation(data_2d[:, 1]) +
                   espressomd.analyze.autocorrelation(data_2d[:, 2]))
        np.testing.assert_allclose(acf, acf_ref, atol=1e-14, rtol=0.)


if __name__ == "__main__":
    ut.main()
