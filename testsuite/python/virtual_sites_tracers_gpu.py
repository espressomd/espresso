#
# Copyright (C) 2013-2018 The ESPResSo project
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
import espressomd
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from espressomd.interactions import FeneBond
from espressomd.utils import handle_errors

from tests_common import verify_lj_forces
from numpy import random
from virtual_sites_tracers_common import VirtualSitesTracersCommon


required_features = "VIRTUAL_SITES_INERTIALESS_TRACERS", "CUDA"


@ut.skipIf(
    not espressomd.gpu_available() or not espressomd.has_features(
        required_features),
           "Test requires VIRTUAL_SITES_INERTIALESS_TRACERS and a GPU")
class VirtualSitesTracers(ut.TestCase, VirtualSitesTracersCommon):

    def setUp(self):
        self.LBClass = lb.LBFluidGPU

if __name__ == "__main__":
    ut.main()
