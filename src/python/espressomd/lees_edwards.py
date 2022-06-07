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

    """
    Interface to the :ref:`Lees-Edwards boundary conditions`.

    Attributes
    ----------
    protocol : :obj:`object`
        Lees--Edwards protocol.
    shear_velocity: :obj:`float`
        Current shear velocity.
    pos_offset : :obj:`float`
        Current position offset
    shear_direction : :obj:`int`
        Shear direction: 0 for the *x*-axis, 1 for the *y*-axis,
        2 for the *z*-axis.
    shear_plane_normal : :obj:`int`
        Shear plane normal: 0 for the *x*-axis, 1 for the *y*-axis,
        2 for the *z*-axis.

    """

    _so_name = "LeesEdwards::LeesEdwards"


@script_interface_register
class Off(ScriptInterfaceHelper):

    """Lees--Edwards protocol resulting in un-shifted boundaries."""
    _so_name = "LeesEdwards::Off"


@script_interface_register
class LinearShear(ScriptInterfaceHelper):

    """Lees--Edwards protocol for linear shear.

    Parameters
    ----------
    initial_pos_offset : :obj:`float`
       Positional offset at the Lees--Edwards boundary at t=0.
    shear_velocity : :obj:`float`
       Shear velocity (velocity jump) across the Lees--Edwards boundary.

    """
    _so_name = "LeesEdwards::LinearShear"


@script_interface_register
class OscillatoryShear(ScriptInterfaceHelper):

    """Lees--Edwards protocol for oscillatory shear.

    Parameters
    ----------
    initial_pos_offset : :obj:`float`
       Positional offset at the Lees--Edwards boundary at t=0.
    amplitude : :obj:`float`
       Maximum amplitude of the positional offset at the Lees--Edwards boundary.
    omega : :obj:`float`
       Radian frequency of the oscillation.
    time_0 : :obj:`float`
       Time offset of the oscillation.

    """
    _so_name = "LeesEdwards::OscillatoryShear"
