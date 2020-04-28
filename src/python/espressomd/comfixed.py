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


@script_interface_register
class ComFixed(ScriptInterfaceHelper):

    """Fix the center of mass of specific types.

    Subtracts mass-weighted fraction of the total
    force action on all particles of the type from
    the particles after each force calculation. This
    keeps the center of mass of the type fixed iff
    the total momentum of the type is zero.

    Parameters
    ----------
    types : array_like
        List of types for which the center of mass should be fixed.
    """

    _so_name = "ComFixed"
    _so_creation_policy = "GLOBAL"
