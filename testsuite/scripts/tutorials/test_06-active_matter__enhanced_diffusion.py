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
import os

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/06-active_matter/solutions/enhanced_diffusion.py",
    cmd_arguments=[5.0], SAMP_STEPS=500, SAMP_LENGTH=100)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_file_generation(self):
        for name in ["msd_{}_0.dat", "vacf_{}_0.dat", "avacf_{}_0.dat"]:
            filepath = os.path.join(tutorial.outdir, name.format(tutorial.vel))
            self.assertTrue(
                os.path.isfile(filepath),
                filepath + " not created")


if __name__ == "__main__":
    ut.main()
