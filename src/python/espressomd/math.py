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
class CylindricalTransformationParameters(ScriptInterfaceHelper):
    """
    Class to hold and validate the parameters needed for a cylindrical transformation.
    The three parameters are available as attributes but are read-only.

    Parameters
    ----------
    center : (3,) array_like of :obj:`float`, default = [0, 0, 0]
        Position of the origin of the cylindrical coordinate system.
    axis : (3,) array_like of :obj:`float`, default = [0, 0, 1]
        Orientation vector of the ``z``-axis of the cylindrical coordinate system.
    orientation: (3,) array_like of :obj:`float`, default = [1, 0, 0]
        The axis on which ``phi = 0``.

    Notes
    -----
    If you provide no arguments, the defaults above are set.
    If you provide only a ``center`` and an ``axis``, an ``orientation`` will be automatically generated that is orthogonal to ``axis``.
    """
    _so_name = "CylindricalTransformationParameters"
