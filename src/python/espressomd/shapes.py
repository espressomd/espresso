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

import collections.abc

from .script_interface import ScriptInterfaceHelper, script_interface_register, ScriptObjectRegistry
from .utils import requires_experimental_features


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
@requires_experimental_features("No test coverage")
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
    radius: :obj:`float`
        The radius of the pore.
    length: :obj:`float`
        The distance between the planes.
    smoothing_radius: :obj:`float`
        Radius of the torus segments
    axis: (3,) array_like of :obj:`float`
        Axis of the cylinder and normal of the planes
    center: (3,) array_like of :obj:`float`
        Position of the center of the cylinder.

    """
    _so_name = "Shapes::SimplePore"


@script_interface_register
class HollowConicalFrustum(Shape, ScriptInterfaceHelper):
    """
    Hollow conical frustum shape.

    Attributes
    ----------
    r1: :obj:`float`
        Radius r1.
    r2: :obj:`float`
        Radius r2.
    length: :obj:`float`
        Length of the conical frustum along ``axis``.
    axis: (3,) array_like of :obj:`float`
        Symmetry axis.
    center: (3,) array_like of :obj:`float`
        Position of the center.


    .. image:: figures/conical_frustum.png
    """
    _so_name = "Shapes::HollowConicalFrustum"


@script_interface_register
class Union(Shape, ScriptObjectRegistry):
    """A union of shapes.

    This shape represents a union of shapes where the distance to the union
    is defined by the smallest distance to any shape contained in the union.

    """
    _so_name = "Shapes::Union"

    def add(self, shape):
        """
        Add a shape to the union.

        Parameters
        ----------
        shape : array_like / instance of :class:`espressomd.shapes.Shape`
            Shape instance(s) to be added to the union.

        """

        def _add(self, shape):
            if isinstance(shape, Shape):
                self.call_method("add", shape=shape)
            else:
                raise ValueError("Only shapes can be added.")

        if isinstance(shape, collections.abc.Iterable):
            for s in shape:
                _add(self, s)
        else:
            _add(self, shape)

    def remove(self, shape):
        """
        Remove a shape from the union.

        Parameters
        ----------
        shape : array_like / instance of :class:`espressomd.shapes.Shape`
            Shape instance(s) to be removed from the union.

        """

        def _remove(self, shape):
            if isinstance(shape, Shape):
                self.call_method("remove", shape=shape)
            else:
                raise ValueError("Only shapes can be removed.")
        if isinstance(shape, collections.abc.Iterable):
            for s in shape:
                _remove(self, s)
        else:
            _remove(self, shape)

    def clear(self):
        """
        Remove all shapes from the union.

        """
        self.call_method("clear")

    def size(self):
        """
        Number of shapes contained in the union.

        """
        return self.call_method("size")
