#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd
import espressomd.electrostatics
import unittest as ut
import unittest_decorators as utx
import tests_common


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures("ELECTROSTATICS")
class P3MGPU_test(ut.TestCase):

    def test(self):
        es = espressomd.System(box_l=[10.0, 10.0, 10.0])
        test_params = {}
        test_params["prefactor"] = 2
        test_params["cao"] = 2
        test_params["r_cut"] = 0.9
        test_params["accuracy"] = 1e-1
        test_params["mesh"] = [10, 10, 10]
        test_params["epsilon"] = 20.0
        test_params["alpha"] = 1.1
        test_params["tune"] = False

        p3m = espressomd.electrostatics.P3MGPU(**test_params)
        es.actors.add(p3m)
        tests_common.assert_params_match(self, test_params, p3m.get_params())


if __name__ == "__main__":
    ut.main()
