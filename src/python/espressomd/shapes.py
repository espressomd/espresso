from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class Cylinder(ScriptInterfaceHelper):
    """
    A cylinder shape.

    Attributes
    ----------
    center : array_like :obj:`float`
             Coordinates of the center of the cylinder.
    axis : array_like :obj:`int`
           Axis of the cylinder.
    direction : :obj:`int`
                Surface orientation, for +1 the normal points
                out of the mantel, for -1 it points inward.
    radius : :obj:`float`
             Radius of the cylinder.
    length : :obj:`float`
             Length of the cylinder.

    """
    _so_name = "Shapes::Cylinder"

@script_interface_register
class HollowCone(ScriptInterfaceHelper):
    """
    A hollow cone shape.

    Attributes
    ----------
    inner_radius : :obj:`float`
                   Inner radius of the cone.
    outer_radius  : :obj:`float`
                    Outer radius of the cone.
    opening_angle  : :obj:`float`
                     Opening angle of the cone (in rad).
    orientation_x  : :obj:`float`
                     x component of the orientation of the cone.
    orientation_y  : :obj:`float`
                     y component of the orientation of the cone.
    orientation_z  : :obj:`float`
                     z component of the orientation of the cone.
    position_x  : :obj:`float`
                  x component of the position of the cone.
    position_y  : :obj:`float`
                  y component of the position of the cone.
    position_z  : :obj:`float`
                  z component of the position of the cone.
    width : :obj:`float`
            Wall thickness of the cone.
    direction : :obj:`int`
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
    cylrad : :obj:`float`
             Radius of the tubes.
    dim : :obj:`float`
          Dimension of the maze.
    nsphere : :obj:`int`
              Number of spherical cavities.
    sphrad : :obj:`float`
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
    axis : array_like :obj:`int`
           Orientation of the pore.
    length : :obj:`float`
             Length of the pore.
    outer_rad_left : :obj:`float`
                     Radius of the left (with respect to the axis) bounding
                     cylinder.
    outer_rad_right : :obj:`float`
                      Radius of the right (with respect to the axis) bounding
                      cylinder.
    pos : array_like :obj:`float`
          Position of the center of the pore
    rad_left : :obj:`float`
               Radius of the left (with respect to the axis) opening.
    rad_right : :obj:`float`
                Radius of the right (with respect to the axis) opening.
    smoothing_radius : :obj:`float`
                       Radius of the smoothing at the opening.

    """
    _so_name = "Shapes::Pore"


@script_interface_register
class Rhomboid(ScriptInterfaceHelper):
    """
    An parallelepiped.

    Attributes
    ----------
    a : array_like :obj:`float`
        First base vector.
    b : array_like :obj:`float`
        Second base vector.
    c : array_like :obj:`float`
        Third base vector.
    corner : array_like :obj:`float`
             Lower left corner of the rhomboid.
    direction : :obj:`int`
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

    channel_width : :obj:`float`
    lower_smoothing_radius : :obj:`float`
    pore_length : :obj:`float`
    pore_mouth : :obj:`float`
    pore_width : :obj:`float`
    upper_smoothing_radius : :obj:`float`

    """
    _so_name = "Shapes::Slitpore"


@script_interface_register
class Sphere(ScriptInterfaceHelper):
    """
    A sphere.

    Attributes
    ----------
    center : array_like :obj:`float`
             Center of the sphere
    radius : :obj:`float`
             Radius of the sphere.
    direction : :obj:`int`
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
    center : array_like :obj:`float`
             Coordinates of the center of the cylinder.
    axis : array_like :obj:`int`
           Axis of the cylinder.
    direction : :obj:`int`
                Surface orientation, for +1 the normal points
                out of the mantel, for -1 it points inward.
    radius : :obj:`float`
             Radius of the cylinder.
    length : :obj:`float`
             Length of the cylinder (not including the caps).

    """
    _so_name = "Shapes::SpheroCylinder"

@script_interface_register
class Stomatocyte(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    inner_radius : :obj:`float`
                   Inner radius of the cone.
    outer_radius : :obj:`float`
                   Outer radius of the cone.
    orientation_x : :obj:`float`
                    x component of the orientation of the cone.
    orientation_y : :obj:`float`
                    y component of the orientation of the cone.
    orientation_z : :obj:`float`
                    z component of the orientation of the cone.
    position_x : :obj:`float`
                 x component of the position of the cone.
    position_y : :obj:`float`
                 y component of the position of the cone.
    position_z : :obj:`float`
                 z component of the position of the cone.
    layer_width : :obj:`float`
    direction : :obj:`int`
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
    dist : :obj:`float`
           Distance from the origin.
    normal : array_like :obj:`int`
             Normal vector of the plan (needs not to be length 1).

    """
    _so_name = "Shapes::Wall"

@script_interface_register
class SimplePore(ScriptInterfaceHelper):
    """
    Two parallel infinite planes, and a cylindrical orfice connecting them.
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
    axis: array_like
       Axis of the cylinder and normal of the planes
    center: array_like
       Position of the center of the cylinder.

    """
    _so_name = "Shapes::SimplePore"
