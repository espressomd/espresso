#
# Copyright (C) 2021-2026 The ESPResSo project
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
    When writing H5MD files, the shear direction and shear plane normals
    are written as integers instead of characters, with 0 = *x*-axis,
    1 = *y*-axis, 2 = *z*-axis.

    Attributes
    ----------
    protocol : :obj:`object`
        Lees--Edwards protocol.
    shear_velocity : :obj:`float`, read-only
        Current shear velocity (derived from the protocol and simulation time).
    pos_offset : :obj:`float`, read-only
        Current position offset (derived from the protocol and simulation time).
    shear_direction : :obj:`str`, {'x', 'y', 'z'}, read-only
        Shear direction. Set via :meth:`set_boundary_conditions`.
    shear_plane_normal : :obj:`str`, {'x', 'y', 'z'}, read-only
        Shear plane normal. Set via :meth:`set_boundary_conditions`.

    Methods
    -------
    set_boundary_conditions()
        Set a protocol, the shear direction and shear normal.
        Must be called at least once before assigning ``protocol``
        directly. Pass ``protocol=None`` to fully disable
        Lees--Edwards (in which case ``shear_direction`` and
        ``shear_plane_normal`` are ignored).

        Parameters
        ----------
        protocol : :obj:`object` or ``None``
        shear_direction : :obj:`str`, {'x', 'y', 'z'}
        shear_plane_normal : :obj:`str`, {'x', 'y', 'z'}

        Raises
        ------
        ValueError
            If ``shear_direction`` or ``shear_plane_normal`` is invalid.
        ValueError
            If ``shear_direction == shear_plane_normal``.

    """

    _so_name = "LeesEdwards::LeesEdwards"
    _so_bind_methods = ("set_boundary_conditions",)


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
       Positional offset at the Lees--Edwards boundary at ``time=time_0``.
    shear_velocity : :obj:`float`
       Shear velocity (velocity jump) across the Lees--Edwards boundary.
    time_0 : :obj:`float`
       Reference time for the linear shear. The positional offset is
       ``initial_pos_offset + (time - time_0) * shear_velocity``.

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
    decay_rate : :obj:`float`
       Apply an exponential decay with the given rate to the
       oscillation's amplitude. Defaults to 0, i.e., no decay.

    """
    _so_name = "LeesEdwards::OscillatoryShear"
