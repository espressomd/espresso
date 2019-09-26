# Copyright (C) 2010-2019 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .__init__ import has_features
import espressomd.lbboundaries


if any(has_features(i) for i in ["LB_BOUNDARIES", "LB_BOUNDARIES_GPU"]):
    @script_interface_register
    class EKBoundaries(espressomd.lbboundaries.LBBoundaries):

        """
        Creates a set of electrokinetics boundaries.

        """
        pass

    @script_interface_register
    class EKBoundary(espressomd.lbboundaries.LBBoundary):

        """
        Creates a EK boundary.

        """
        pass
