#
# Copyright (C) 2019-2020 The ESPResSo project
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
import unittest_decorators as utx
import stokesian_dynamics as sd


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsSetupTest(sd.StokesianDynamicsSetupTest):
    device = 'cpu'

    def test_pbc_checks(self):
        self.pbc_checks()


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsTest(sd.StokesianDynamicsTest):
    device = 'cpu'

    def test_default(self):
        self.falling_spheres(1.0, 1.0, 1.0)

    def test_rescaled(self):
        self.falling_spheres(1.0, 4.5, 2.5)

    def test_different_time_step(self):
        self.falling_spheres(0.7, 1.0, 1.0)


@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDiffusionTest(sd.StokesianDiffusionTest):
    device = 'cpu'

    def test(self):
        self.check()


if __name__ == '__main__':
    ut.main()
