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
import espressomd.observables


def disable_GUI(code):
    # integrate without visualizer
    breakpoint = "visualizer.run(1)"
    assert breakpoint in code
    code = code.replace(breakpoint, "steps=1\nsystem.integrator.run(steps)", 1)
    return code


sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "/home/thilo/code/espresso2/espresso/testsuite/scripts/samples/local_samples/visualization_poiseuille.py",
    substitutions=disable_GUI, steps=100)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_particle_coupling(self):
        part_vel = espressomd.observables.ParticleVelocities(
            ids=list(range(100)))
        mean_velocity = np.mean(part_vel.calculate())
        self.assertGreater(mean_velocity, 1e-5)


if __name__ == "__main__":
    ut.main()
