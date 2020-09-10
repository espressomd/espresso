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

import unittest as ut
import importlib_wrapper
import numpy as np

try:
    import pint  # pylint: disable=unused-import
except ImportError:
    tutorial = importlib_wrapper.MagicMock()
    skipIfMissingFeatures = ut.skip(
        "Python module pint not available, skipping test!")
else:
    tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
        "@TUTORIALS_DIR@/13-charged_rod/13-charged_rod.py")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_equilibration(self):
        # TODO 
        # this is just a placeholder
        self.assertLess(tutorial.energies[-1], tutorial.energies[0])

    def test_file_generation(self):
        # test .vsf/.vcf files exist
        for name in ["ion_condensation.vcf", "ion_condensation.vsf"]:
            filepath = os.path.join(tutorial.outdir, name)
            self.assertTrue(
                os.path.isfile(filepath),
                filepath + " not created")


if __name__ == "__main__":
    ut.main()
