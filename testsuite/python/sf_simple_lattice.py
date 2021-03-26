#
# Copyright (C) 2017-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published byss
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
import espressomd
import numpy as np
from itertools import product


class sfTest(ut.TestCase):

    box_l = 20
    part_ty = 0
    sf_order = 20
    s = espressomd.System(box_l=[box_l, box_l, box_l])

    def make_lattice(self):
        xen, yen, zen = [
            x for x in range(0, self.box_l, 1)], [
            x for x in range(0, self.box_l, 2)], [
            x for x in range(0, self.box_l, 4)]
        for i, j, k in product(xen, yen, zen):
            self.s.part.add(type=self.part_ty, pos=(i, j, k))

    def peaks(self):
        a, b, c = [0, 2 * np.pi], [0, np.pi], [0, np.pi / 2]
        self.diags = [np.sqrt(a**2 + b**2 + c**2)
                      for (a, b, c) in product(a, b, c)]
        self.diags.remove(0)

    def scatter_lattice(self):
        sf_observable_x, sf_observable_y = self.s.analysis.structure_factor(
            sf_types=[self.part_ty], sf_order=self.sf_order)
        self.sf_data = [(x, y)
                        for (x, y) in zip(sf_observable_x, sf_observable_y)]
        self.sf_data.sort(key=lambda tup: tup[1])

    def test(self):
        self.make_lattice()
        self.peaks()
        self.scatter_lattice()
        self.assertTrue([np.any(np.isclose(element, self.sf_data, rtol=1e-02))
                         for element in self.diags])


if __name__ == "__main__":
    ut.main()
