#
# Copyright (C) 2010-2022 The ESPResSo project
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

import math
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
    If you provide only a ``center`` and an ``axis``, an ``orientation``
    will be automatically generated that is orthogonal to ``axis``.

    """
    _so_name = "Math::CylindricalTransformationParameters"


def calc_quaternions_from_angles(polar_angle, azimuthal_angle):
    """
    Convert the spherical coordinate representation of a vector
    to a quaternion that corresponds to that direction.
    Note that this quaternion is not unique, as the third degree of freedom
    (rotation around the given vector) is not specified.

    Parameters
    ----------
    polar_angle : :obj:`float`
        The polar angle (angle with the z-axis) in the range [0, pi].
    azimuthal_angle : :obj:`float`
        The azimuthal angle (angle with the x-axis) in the range [0, 2 pi).

    Returns
    -------
    quaternion : (4,) list of :obj:`float`

    """

    q1w = math.cos(polar_angle / 2.)
    q1x = 0.
    q1y = math.sin(polar_angle / 2.)
    q1z = 0.

    q2w = math.cos(azimuthal_angle / 2.)
    q2x = 0.
    q2y = 0.
    q2z = math.sin(azimuthal_angle / 2.)

    q3w = q1w * q2w - q1x * q2x - q1y * q2y - q1z * q2z
    q3x = q1w * q2x + q1x * q2w - q1y * q2z + q1z * q2y
    q3y = q1w * q2y + q1x * q2z + q1y * q2w - q1z * q2x
    q3z = q1w * q2z - q1x * q2y + q1y * q2x + q1z * q2w

    return [q3w, q3x, q3y, q3z]
