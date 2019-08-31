# Copyright (C) 2010-2018 The ESPResSo project
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


class Shape:
    _so_bind_methods = ("calc_distance",)


@script_interface_register
class Cylinder(Shape, ScriptInterfaceHelper):

    """
    A cylinder shape.

    Attributes
    ----------
    center : (3,) array_like of :obj:`float`
        Coordinates of the center of the cylinder.
    axis : (3,) array_like of :obj:`float`
        Axis of the cylinder.
    radius : :obj:`float`
        Radius of the cylinder.
    length : :obj:`float`
        Length of the cylinder.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.
    open : :obj:`bool`
        cylinder is open or has caps.

    """
    _so_name = "Shapes::Cylinder"


@script_interface_register
class Ellipsoid(Shape, ScriptInterfaceHelper):

    """
    An ellipsoid.

    For now only ellipsoids of revolution are supported.
    The symmetry axis is aligned parallel to the x-direction.

    Attributes
    ----------
    center : (3,) array_like of :obj:`float`
        Coordinates of the center of the ellipsoid.
    a : :obj:`float`
        Semiaxis along the axis of rotational symmetry.
    b : :obj:`float`
        Equatorial semiaxes.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.
    """
    _so_name = "Shapes::Ellipsoid"


@script_interface_register
class HollowCone(Shape, ScriptInterfaceHelper):

    """
    A hollow cone shape.

    Attributes
    ----------
    inner_radius : :obj:`float`
        Inner radius of the cone.
    outer_radius  : :obj:`float`
        Outer radius of the cone.
    opening_angle : :obj:`float`
        Opening angle of the cone (in rad).
    axis : (3,) array_like of :obj:`float`
        Axis of symmetry, prescribes orientation of the cone.
    center : (3,) array_like of :obj:`float`
        Position of the cone.
    width : :obj:`float`
        Wall thickness of the cone.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.

    """
    _so_name = "Shapes::HollowCone"


@script_interface_register
class Rhomboid(Shape, ScriptInterfaceHelper):

    """
    An parallelepiped.

    Attributes
    ----------
    a : (3,) array_like of :obj:`float`
        First base vector.
    b : (3,) array_like of :obj:`float`
        Second base vector.
    c : (3,) array_like of :obj:`float`
        Third base vector.
    corner : (3,) array_like of :obj:`float`
        Lower left corner of the rhomboid.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.

    """
    _so_name = "Shapes::Rhomboid"


@script_interface_register
class Slitpore(Shape, ScriptInterfaceHelper):

    """

    .. image:: figures/slitpore.png

    Attributes
    ----------

    channel_width : :obj:`float`
    lower_smoothing_radius : :obj:`float`
    pore_length : :obj:`float`
    pore_mouth : :obj:`float`
    pore_width : :obj:`float`
    upper_smoothing_radius : :obj:`float`
    dividing_plane : :obj:`float`

    """
    _so_name = "Shapes::Slitpore"


@script_interface_register
class Sphere(Shape, ScriptInterfaceHelper):

    """
    A sphere.

    Attributes
    ----------
    center : (3,) array_like of :obj:`float`
        Center of the sphere
    radius : :obj:`float`
        Radius of the sphere.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.

    """
    _so_name = "Shapes::Sphere"


@script_interface_register
class SpheroCylinder(Shape, ScriptInterfaceHelper):

    """
    A cylinder with hemispheres as caps.

    Attributes
    ----------
    center : (3,) array_like of :obj:`float`
        Coordinates of the center of the cylinder.
    axis : (3,) array_like of :obj:`float`
        Axis of the cylinder.
    radius : :obj:`float`
        Radius of the cylinder.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.
    length : :obj:`float`
        Length of the cylinder (not including the caps).

    """
    _so_name = "Shapes::SpheroCylinder"


@script_interface_register
class Stomatocyte(Shape, ScriptInterfaceHelper):

    """
    Attributes
    ----------
    inner_radius : :obj:`float`
        Inner radius of the stomatocyte.
    outer_radius : :obj:`float`
        Outer radius of the stomatocyte.
    axis : (3,) array_like of :obj:`float`
        Symmetry axis, prescribing the orientation of the stomatocyte.
    center : (3,) array_like of :obj:`float`
        Position of the stomatocyte.
    layer_width : :obj:`float`
        Scaling parameter.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.

    """

    _so_name = "Shapes::Stomatocyte"


@script_interface_register
class Torus(Shape, ScriptInterfaceHelper):

    """
    A torus shape.
    Attributes
    ----------
    center : (3,) array_like of :obj:`float`
        Coordinates of the center of the torus.
    normal : (3,) array_like of :obj:`float`
        Normal axis of the torus.
    radius : :obj:`float`
        Radius of the torus.
    tube_radius : :obj:`float`
        Radius of the tube.
    direction : :obj:`int`
        Surface orientation, for +1 the normal points
        out of the mantel, for -1 it points inward.

    """
    _so_name = "Shapes::Torus"


@script_interface_register
class Wall(Shape, ScriptInterfaceHelper):

    """
    An infinite plane.

    Attributes
    ----------
    dist : :obj:`float`
        Distance from the origin.
    normal : (3,) array_like of :obj:`int`
        Normal vector of the plane (needs not to be length 1).

    """
    _so_name = "Shapes::Wall"


@script_interface_register
class SimplePore(Shape, ScriptInterfaceHelper):

    """
    Two parallel infinite planes, and a cylindrical channel connecting them.
    The cylinder and the planes are connected by torus segments with an
    adjustable radius.

    Attributes
    ----------
    radius: float
        The radius of the pore.
    length: float
        The distance between the planes.
    smoothing_radius: float
        Radius of the torus segments
    axis: (3,) array_like of :obj:`float`
        Axis of the cylinder and normal of the planes
    center: (3,) array_like of :obj:`float`
        Position of the center of the cylinder.

    """
    _so_name = "Shapes::SimplePore"
