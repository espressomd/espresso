#
# Copyright (C) 2025-2026 The ESPResSo project
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
import importlib_wrapper
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/boltzmann_inversion/boltzmann_inversion.py",
    p3m_params={"mesh": [16, 16, 16], "cao": 3})


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_boltzmann_inversion(self):
        rdf_exp = tutorial.rdf
        rdf_imp = tutorial.rdf_dh
        np.testing.assert_allclose(rdf_exp, rdf_imp, rtol=0.1, atol=0.25)


if __name__ == "__main__":
    ut.main()
