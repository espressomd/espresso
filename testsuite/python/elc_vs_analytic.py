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
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.electrostatics
from espressomd import electrostatic_extensions


@utx.skipIfMissingFeatures(["P3M"])
class ELC_vs_analytic(ut.TestCase):
    # Handle to espresso system
    box_l = 200.
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    accuracy = 1e-7
    check_accuracy = 1e-4
    elc_gap = 75.0
    system.time_step = 0.01
    delta_mid_top = 1.
    delta_mid_bot = 39. / 41.
    distance = 1.

    number_samples = 25
    zPos = np.linspace(0.1, 8, number_samples)
    q = np.arange(-5.0, 5.1, 2.5)

    def test_elc(self):
        """
        Testing ELC against the analytic solution for an infinite large simulation box with dielectric contrast on the bottom of the box, which can be calculated analytically with image charges.
        """
        self.system.part.add(id=1, pos=self.system.box_l / 2., q=self.q[0])
        self.system.part.add(id=2, pos=self.system.box_l / 2. + [0, 0, 1],
                             q=-self.q[0])

        self.system.box_l = [self.box_l, self.box_l, self.box_l + self.elc_gap]
        self.system.cell_system.set_domain_decomposition(
            use_verlet_lists=True)
        self.system.periodicity = [1, 1, 1]
        p3m = espressomd.electrostatics.P3M(prefactor=1.,
                                            accuracy=self.accuracy,
                                            mesh=[58, 58, 70],
                                            cao=4)
        self.system.actors.add(p3m)

        elc = electrostatic_extensions.ELC(gap_size=self.elc_gap,
                                           maxPWerror=self.accuracy,
                                           delta_mid_bot=self.delta_mid_bot,
                                           delta_mid_top=self.delta_mid_top)
        self.system.actors.add(elc)

        elc_results = self.scan()

        # ANALYTIC SOLUTION
        charge_reshaped = np.square(self.q.reshape(-1, 1))
        analytic_force = charge_reshaped * (1 / self.distance ** 2 + self.delta_mid_bot * (
            1 / np.square(2 * self.zPos) - 1 / np.square(2 * self.zPos + self.distance)))
        analytic_energy = charge_reshaped * (-1 / self.distance + self.delta_mid_bot * (1 / (
            4 * self.zPos) - 1 / (2 * self.zPos + self.distance) + 1 / (4 * (self.zPos + self.distance))))

        analytic_results = np.dstack((analytic_force, analytic_energy))

        np.testing.assert_allclose(
            elc_results, analytic_results, rtol=0, atol=self.check_accuracy)

    def scan(self):
        result_array = np.empty((len(self.q), len(self.zPos), 2))
        for chargeIndex, charge in enumerate(self.q):
            self.system.part[1].q = charge
            self.system.part[2].q = -charge
            for i, z in enumerate(self.zPos):
                pos = self.system.part[1].pos
                self.system.part[1].pos = [pos[0], pos[1], z]
                self.system.part[2].pos = [pos[0], pos[1], z + self.distance]

                self.system.integrator.run(0)
                result_array[chargeIndex, i, 0] = self.system.part[1].f[2]
                result_array[chargeIndex, i, 1] = self.system.analysis.energy()[
                    "total"]
        return result_array


if __name__ == "__main__":
    ut.main()
