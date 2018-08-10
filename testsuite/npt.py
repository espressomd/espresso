#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
from espressomd import thermostat
import tests_common


@ut.skipIf(not espressomd.has_features(["NPT", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class NPTintegrator(ut.TestCase):
    """This compares pressure and compressibility of a LJ system against expected values."""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    p_ext = 2.0

    def setUp(self):
        box_l = 5.78072780098
        self.system.box_l = [box_l] * 3
        self.system.time_step = 0.02
        self.system.cell_system.skin = 0.4

        data = np.genfromtxt(tests_common.abspath(
            "data/npt_lj_system.data"))

        # Input format: id pos f
        for particle in data:
            pos = particle[:3]
            f = particle[3:]
            self.system.part.add(pos=pos, f=f)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1,
            cutoff=1.12246, shift=0.25)

        self.system.thermostat.set_npt(kT=1.0, gamma0=1, gammav=0.002)
        self.system.integrator.set_isotropic_npt(
            ext_pressure=self.p_ext, piston=0.001)

    def test_npt(self):
        avp = 0
        n = 600
        Vs = []
        self.system.integrator.run(100)
        for t in range(n):
            self.system.integrator.run(15)
            avp += self.system.analysis.pressure()['total']
            l0 = self.system.box_l[0]
            Vs.append(pow(l0, 3))

        avp /= n
        compressibility = pow(np.std(Vs), 2) / np.average(Vs)

        self.assertAlmostEqual(2.0, avp, delta=0.02)
        self.assertAlmostEqual(0.2, compressibility, delta=0.013)


if __name__ == "__main__":
    ut.main()
