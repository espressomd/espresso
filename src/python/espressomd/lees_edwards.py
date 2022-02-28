#
# Copyright (C) 2021-2022 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class LeesEdwards(ScriptInterfaceHelper):

    """Interface to the Lees-Edwards boundary conditions.

       See documentation.

    """

    _so_name = "LeesEdwards::LeesEdwards"


@script_interface_register
class Off(ScriptInterfaceHelper):

    """Lees-Edwards protocol resulting in un-shifted boundaries."""
    _so_name = "LeesEdwards::Off"


@script_interface_register
class LinearShear(ScriptInterfaceHelper):

    """Lees-Edwards protocol for linear shear.

    Parameters
    ----------
    shear_direction
       Cartesian coordinate of the shear direction (0=x, 1=y, 2=z)
    shear_plane_normal
       Cartesian coordinate of the shear plane normal
    initial_pos_offset
       Positional offset at the Lees-Edwards boundary at t=0
    shear_velocity
       Shear velocity (velocity jump) across the Lees-Edwards boundary

    """
    _so_name = "LeesEdwards::LinearShear"


@script_interface_register
class OscillatoryShear(ScriptInterfaceHelper):

    """Lees-Edwards protocol for oscillatory shear.

    Parameters
    ----------
    shear_direction
       Cartesian coordinate of the shear direction (0=x, 1=y, 2=z)
    shear_plane_normal
       Cartesian coordinate of the shear plane normal
    amplitude
       Maximum amplitude of the positional offset at the Lees-Edwards boundary
    frequency
       Frequency of the shear
    time_0
       Time offset of the oscillation


    """
    _so_name = "LeesEdwards::OscillatoryShear"
