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

try:
    import pint  # pylint: disable=unused-import
except ImportError:
    tutorial = importlib_wrapper.MagicMock()
    skipIfMissingFeatures = ut.skip(
        "Python module pint not available, skipping test!")
else:
    tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
        "/home/thilo/code/espresso2/espresso/testsuite/scripts/tutorials/local_tutorials/constant_pH/constant_pH.py", script_suffix="interactions",
        USE_WCA=True, USE_ELECTROSTATICS=True, NUM_PHS=8, NUM_SAMPLES=10)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test(self):
        pass


if __name__ == "__main__":
    ut.main()
