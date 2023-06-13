# Copyright (C) 2019-2022 The ESPResSo project
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

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "/home/thilo/code/espresso2/espresso/testsuite/scripts/tutorials/local_tutorials/langevin_dynamics/langevin_dynamics.py")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_diffusion_coefficient(self):
        D_gk = tutorial.diffusion_gk
        D_msd = tutorial.diffusion_msd
        D_ref = tutorial.KT / np.array(tutorial.gammas)
        np.testing.assert_allclose(D_msd, D_ref, rtol=0, atol=0.02)
        np.testing.assert_allclose(D_gk, D_ref, rtol=0, atol=0.02)


if __name__ == "__main__":
    ut.main()
