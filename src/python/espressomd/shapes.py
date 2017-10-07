from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
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

@script_interface_register
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

@script_interface_register
class Maze(ScriptInterfaceHelper):
    """
    Spherical cavities on a regular grid that are
    connected by tubes.

    Attributes
    ----------
    cylrad : float
       Radius of the tubes.
    dim : float
       Dimension of the maze.
    nsphere : int
       Number of spherical cavities.
    sphrad : float
       Radius of the spherical cavities.
    """
    _so_name = "Shapes::Maze"

@script_interface_register
class Pore(ScriptInterfaceHelper):
    """
    A cylinder with a conical pore between the faces. The pore openings
    are smoothed with torus segment. The outer radius can be chosen such
    that it is bigger than the box, to get a wall with a pore.

    Attributes
    ----------
    axis : array_like
       Orientation of the pore.
    length : float
       Length of the pore.
    outer_rad_left : float
       Radius of the left (with respect to the axis) bounding cylinder.
    outer_rad_right : float
       Radius of the right (with respect to the axis) bounding cylinder.
    pos : array_like
       Position of the center of the pore
    rad_left : float
       Radius of the left (with respect to the axis) opening.
    rad_right : float
       Radius of the right (with respect to the axis) opening.
    smoothing_radius : float
       Radius of the smoothing at the opening.
    """
    _so_name = "Shapes::Pore"


@script_interface_register
class Rhomboid(ScriptInterfaceHelper):
    """
    An parallelepiped.

    Attributes
    ----------
    a : array_like
       First base vector.
    b : array_like
       Second base vector.
    c : array_like
       Third base vector.
    corner : array_like
       Lower left corner of the rhomboid.
    direction : int
       Surface orientation, for +1 the normal points
       out of the mantel, for -1 it points inward.
    """
    _so_name = "Shapes::Rhomboid"


@script_interface_register
class Slitpore(ScriptInterfaceHelper):
    """

    .. image:: figures/slitpore.png

    Attributes
    ----------

    channel_width : float
    lower_smoothing_radius : float
    pore_length : float
    pore_mouth : float
    pore_width : float
    upper_smoothing_radius : float
    """
    _so_name = "Shapes::Slitpore"


@script_interface_register
class Sphere(ScriptInterfaceHelper):
    """
    A sphere.

    Attributes
    ----------
    center : array_like
       Center of the sphere
    radius : float
       Radius of the sphere.
    direction : int
       Surface orientation, for +1 the normal points
       out of the mantel, for -1 it points inward.
    """
    _so_name = "Shapes::Sphere"


@script_interface_register
class SpheroCylinder(ScriptInterfaceHelper):
    """
    A cylinder with hemispheres as caps.

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
      Length of the cylinder (not including the caps).
    """
    _so_name = "Shapes::SpheroCylinder"

@script_interface_register
class Stomatocyte(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    inner_radius : float
       Inner radius of the cone.
    outer_radius  : float
       Outer radius of the cone.

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

    layer_width : float

    direction : int
       Surface orientation, for +1 the normal points
       out of the mantel, for -1 it points inward.
    """

    _so_name = "Shapes::Stomatocyte"


@script_interface_register
class Wall(ScriptInterfaceHelper):
    """
    An infinite plane.

    Attributes
    ----------
    dist : float
       Distance from the origin.
    normal : array_like
       Normal vector of the plan (needs not to be length 1).
    """
    _so_name = "Shapes::Wall"
