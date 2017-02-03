from .script_interface import ScriptInterfaceHelper

class Cylinder(ScriptInterfaceHelper):
    """
    A cylinder shape.

    Attributes
    ----------
    center : array_like
       Coordinates of the center of the cylinder

    axis : array_like
       Axis of the cylinder

    direction : int
       Surface orientation, for +1 the normal points
       out of the mantel, for -1 it points inward.

    radius : float
       Radius of the cylinder.

    length : float
      Length of the cylinder.
    """
    _so_name = "Shapes::Cylinder"

class HollowCone(ScriptInterfaceHelper):
    """
    A hollow cone shape.

    Attributes
    ----------
    inner_radius : float
       Inner radius of the cone.
    outer_radius  : float
       Outer radius of the cone.

    opening_angle  : float
       Opening angle of the cone (in rad)

    orientation_x  : float
       x component of the orientation of the cone.
    orientation_y  : float
       y component of the orientation of the cone.
    orientation_z  : float
       z component of the orientation of the cone.

    position_x  : float
       x component of the position of the cone.
    position_y  : float
       y component of the position of the cone.
    position_z  : float
       z component of the position of the cone.

    width : float
       wall thickness of the cone.

    direction : int
       Surface orientation, for +1 the normal points
       out of the mantel, for -1 it points inward.
    """
    _so_name = "Shapes::HollowCone"

class Maze(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    cylrad : float
    dim : float
    nsphere : float
    sphrad : float
    """
    _so_name = "Shapes::Maze"


class Pore(ScriptInterfaceHelper):
    _so_name = "Shapes::Pore"


class Rhomboid(ScriptInterfaceHelper):
    _so_name = "Shapes::Rhomboid"


class Slitpore(ScriptInterfaceHelper):
    _so_name = "Shapes::Slitpore"


class Sphere(ScriptInterfaceHelper):
    _so_name = "Shapes::Sphere"


class SpheroCylinder(ScriptInterfaceHelper):
    _so_name = "Shapes::SpheroCylinder"

class Stomatocyte(ScriptInterfaceHelper):
    _so_name = "Shapes::Stomatocyte"


class Wall(ScriptInterfaceHelper):
    _so_name = "Shapes::Wall"
