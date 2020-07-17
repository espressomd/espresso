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


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS_GPU"])
class StokesianDynamicsSetupTest(sd.StokesianDynamicsSetupTest):
    device = 'gpu'

    def test_pbc_checks(self):
        self.pbc_checks()


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS_GPU"])
class StokesianDynamicsTest(sd.StokesianDynamicsTest):
    device = 'gpu'

    def test_default_fts(self):
        self.falling_spheres(1.0, 1.0, 1.0, 'fts', sd_short=True)

    def test_default_ft(self):
        self.falling_spheres(1.0, 1.0, 1.0, 'ft', sd_short=True)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS_GPU"])
class StokesianDiffusionTest(sd.StokesianDiffusionTest):
    device = 'gpu'

    def test(self):
        pass


if __name__ == '__main__':
    ut.main()
