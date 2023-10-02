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
import importlib_wrapper as iw
import numpy as np

tutorial, skipIfMissingFeatures = iw.configure_and_import(
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics.py", RUN_TIME=10, TOTAL_FRAMES=10)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_simulation_fields_finite(self):
        for species in (*tutorial.educt_species, *tutorial.product_species):
            assert np.all(np.isfinite(species[:, :, :].density))
            assert np.all(species[:, :, :].density >= 0)

        assert np.all(np.isfinite(tutorial.lbf[:, :, :].velocity))


if __name__ == "__main__":
    ut.main()
