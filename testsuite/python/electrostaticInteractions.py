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
import unittest as ut
import espressomd
import numpy as np
from espressomd.electrostatics import P3M, DH
from tests_common import *


class ElectrostaticInteractionsTests(ut.TestCase):
    
    def setup(self):
        self.system.box_l = 10, 10, 10
        if not self.system.part.exists(0):
            self.system.part.add(id=0, pos=(0.0, 0.0, 0.0), q=1)
        if not self.system.part.exists(1):
            self.system.part.add(id=1, pos=(0.1, 0.1, 0.1), q=-1)

        
    test_P3M = generate_test_for_class(P3M, dict(bjerrum_length=1.0,
                                                                 epsilon=0.0,
                                                                 inter=1000,
                                                                 mesh_off=[
                                                                     0.5, 0.5, 0.5],
                                                                 r_cut=2.4,
                                                                 mesh=[
                                                                     2, 2, 2],
                                                                 cao=1,
                                                                 alpha=12,
                                                                 accuracy=0.01))

    test_DH = generate_test_for_class(DH, dict(bjerrum_length=1.0,
                                                               kappa=2.3,
                                                               r_cut=2))
    if "CDG" in espressomd.features():
        test_CDH = generate_test_for_class(CDH, dict(bjerrum_length=1.0,
                                                                     kappa=2.3,
                                                                     r_cut=2,
                                                                     r0=1,
                                                                     r1=2,
                                                                     eps_int=0.8,
                                                                     eps_ext=1,
                                                                     alpha=2))


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()

