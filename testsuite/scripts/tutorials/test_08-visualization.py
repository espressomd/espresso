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


def disable_visualizer_GUI(code):
    breakpoint = "t = Thread(target=main)"
    assert breakpoint in code
    code = code.split(breakpoint, 1)[0] + "main()"
    return code


tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/08-visualization/08-visualization.py",
    substitutions=disable_visualizer_GUI, int_n_times=5, int_steps=100,
    matplotlib_notebook=False)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system


if __name__ == "__main__":
    ut.main()
