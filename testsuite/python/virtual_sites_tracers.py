#
# Copyright (C) 2013-2022 The ESPResSo project
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
import espressomd.lb

from virtual_sites_tracers_common import VirtualSitesTracersCommon


@utx.skipIfMissingFeatures(
    ["VIRTUAL_SITES_INERTIALESS_TRACERS", "WALBERLA"])
class VirtualSitesTracers(VirtualSitesTracersCommon, ut.TestCase):

    LBClass = espressomd.lb.LBFluidWalberla


if __name__ == "__main__":
    ut.main()
