.. _Single particle forces (constraints):

Single particle forces (constraints)
====================================

:class:`espressomd.constraints.Constraint`

A Constraint is an immobile surface which can interact with particles via a
nonbonded potential, where the distance between the two particles is
replaced by the distance of the center of the particle to the surface.

The constraints are identified like a particle via its type ``particle_type`` for the
non-bonded interaction. After a type is defined for each constraint one
has to define the interaction of all different particle types with the
constraint using the  :class:`espressomd.interactions.NonBondedInteractions` class.

.. _Shaped-based constraints:

Shaped-based constraints
------------------------

.. note::
    `Feature CONSTRAINTS required`

In order to use shapes you first have to import the :mod:`espressomd.shapes`
module. This module provides classes for the different available shapes::

    import espressomd.shapes

Shapes define geometries which can be used in |es| either as
constraints in particle interactions or as a boundary for a
Lattice-Boltzmann fluid. 

To avoid unexpected behaviour make sure all parts of your shape are 
within the central box since the distance to the shape is calculated only 
within the central box. If parts of the shape are placed 
outside of the central box these parts are truncated by the box boundaries. This may 
or may not be desired as for example in the case of a cylinder without or with cylinder cover. 

A shape is instantiated by calling its constructor. If you wanted to
create a wall shape you could do::

    wall = espressomd.shapes.Wall()

Available shapes are listed below.

    - :class:`espressomd.shapes.Cylinder`
    - :class:`espressomd.shapes.HollowCone`
    - :class:`espressomd.shapes.Maze`
    - :class:`espressomd.shapes.Pore`
    - :class:`espressomd.shapes.Rhomboid`
    - :class:`espressomd.shapes.Slitpore`
    - :class:`espressomd.shapes.Sphere`
    - :class:`espressomd.shapes.SpheroCylinder`
    - :class:`espressomd.shapes.Stomatocyte`


.. _Adding shape-based constraints to the system:

Adding shape-based constraints to the system
--------------------------------------------

Usually you want to use constraints based on a shape.
The module :mod:`espressomd.constraints` provides the class
:class:`espressomd.constraints.ShapeBasedConstraint`::

    shape_constraint = espressomd.constraints.ShapeBasedConstraint(shape=my_shape)

In order to add the constraint to to the system
invoke the :meth:`espressomd.constraints.add` method::

    system.constraints.add(shape_constraint)

All previosly listed shapes can be added to the system's constraints 
by passing a initialized shape object to :meth:`system.constraints.add`, returning a constraint object ::
  
    myShape = Wall( dist=20, normal=[0.1, 0.0, 1] )
    myConstraint = system.constraints.add(shape = myShape, particle_type=p_type)

The extra argument ``particle_type`` specifies the nonbonded interaction to be used with
that constraint.

There are two further optional parameters and that can
be used to fine tune the behavior of the constraint. If ``penetrable`` is
set to ``True`` then particles can move through the constraint in this case the
other option ``only_positive`` controls whether the particle is subject to the interaction
potential of the wall. If set to then the constraint will only act in
the direction of the normal vector.

If we wanted to add a non-penetrable pore constraint to our simulation,
we could do the following::

    pore = espressomd.shapes.Pore(axis=[1,0,0], length=2, pos=[15,15,15], smoothing_radius=0.5)
    pore_constraint = espressomd.constraints.ShapeBasedConstraint(shape=pore, penetrable=0, particle_type=1)
    system.constraints.add(pore_constraint)

Interactions between the pore and other particles are then defined
as usual (:ref:`Non-bonded interactions`).

.. _Deleting a constraint:

Deleting a constraint
~~~~~~~~~~~~~~~~~~~~~

Constraints can be removed in a similar fashion using :meth:`espressomd.system.constraints.remove` ::

    system.constraints.remove(myConstraint)

This command will delete the specified constraint.


.. _Getting the currently defined constraints:

Getting the currently defined constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can interate through constraints, for example ::
  
    >>> for c in system.constraints:
    >>>    print(c.shape)

will print the shape information for all defined constraints.


.. _Getting the force on a constraint:

Getting the force on a constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`espressomd.system.constraints.total_force`

Returns the force acting on the a constraint. Note, however, that this
are only forces due to interactions with particles, not with other
constraints. Also, these forces still do not mean that the constraints
move, they are just the negative of the sum of forces acting on all
particles due to this constraint. Similarly, the total energy does not
contain constraint-constraint contributions.

For example the pressure from wall ::
    >>> p = system.constraints[0].total_force()
    >>> print(p)

.. _Getting the minimal distance to a constraint:

Getting the minimal distance to a constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:todo: `This feature is not yet implemented .`

calculates the smallest distance to all non-penetrable
constraints, that can be repulsive (wall, cylinder, sphere, rhomboid,
maze, pore, slitpore). Negative distances mean that the position is
within the area that particles should not access. Helpful to find
initial configurations.)

.. _Available Shapes:

Available Shapes
~~~~~~~~~~~~~~~~

:class:`espressomd.shapes`

Python Syntax::

    import espressomd from espressomd.shapes import <SHAPE>
    system=espressomd.System()

``<SHAPE>`` can be any of the available shapes.

The surface's geometry is defined via a few available shapes.
The following shapes can be used as constraints.

.. warning::
   When using shapes with concave edges and corners, the fact that a particle
   only interacts with the closest point on the constraint surface leads to discontinuous
   force fields acting on the particles. This breaks energy conservation in otherwise
   symplectic integrators. Often, the total energy of the system increases exponentially.


:class:`espressomd.shapes.Wall`
    An infinite plane`.

The resulting surface is a plane defined by the normal vector ``normal`` 
and the distance ``dist`` from the origin (in the direction of the normal vector).
The force acts in direction of the normal. 
Note that ``dist`` describes the distance from the origin in units of the normal 
vector so that the product of ``dist`` and ``normal`` is a point on the surface.
Therefore negative distances are quite common!

.. figure:: figures/shape-wall.png
   :alt: Example constraint with a ``Wall`` shape.
   :align: center
   :height: 6.00000cm
   
Pictured is an example cosntraint with a ``Wall`` shape created with ::

    wall = Wall( dist=20, normal=[0.1,0.0,1] )
    system.constraints.add(shape=wall, particle_type=0)
    
In variant (1) if the only_positive flag is set to 1, interactions are only calculated if
the particle is on the side of the wall in which the normal vector is
pointing.
This has only an effect for penetrable walls. If the flag is
set to 1, then slip boundary interactions apply that are essential for
microchannel flows like the Plane Poiseuille or Plane Couette Flow.
Youalso need to use the tunable\_slip interaction (see [sec:tunableSlip])
for this too work.


:class:`espressomd.shapes.Sphere`
    A sphere.

The resulting surface is a sphere with center ``center`` and radius ``radius``. 
The direction ``direction`` determines the force direction, ``-1`` or for inward and ``+1`` for outward.

.. _shape-sphere:

.. figure:: figures/shape-sphere.png
   :alt: Example constraint with a ``Sphere`` shape.
   :align: center
   :height: 6.00000cm
   
Pictured is an example cosntraint with a ``Sphere`` shape created with ::
  
    sphere = Sphere(center=[25,25,25], radius = 15, direction = 1 )
    system.constraints.add(shape=sphere, particle_type=0)


:class:`espressomd.shapes.Cylinder`
    A cylinder

The resulting surface is a cylinder with center ``center`` and radius ``radius``.
The ``length`` parameter is **half** of the cylinder length.
The ``axis`` parameter is a vector along the cylinder axis, which is normalized in the program.
The direction ``direction`` determines the force direction, ``-1`` or for inward and ``+1`` for outward.



.. figure:: figures/shape-cylinder.png
   :alt: Example constraint with a ``Cylinder`` shape.
   :align: center
   :height: 6.00000cm
   
Pictured is an example constraint with a ``Cylinder`` shape created with ::

    cylinder=Cylinder(center=[25, 25, 25], axis = [1, 0, 0], direction = 1, radius = 10, length = 30)
    system.constraints.add(shape=cylinder, particle_type = 0)

:class:`espressomd.shapes.Rhomboid`
    A rhomboid or parallelpiped.

:todo: `This shape is currently broken. Please do not use.`

The resulting surface is a rhomboid, defined by one corner located at ``corner`` 
and three adjacent edges, defined by the three vectors connecting the 
corner ``corner`` with it’s three neighboring corners:
``a`` ``[ax ay az ]``; ``b`` ``[bx by bz]`` and ``c`` ``[cx cy cz]``.
The direction ``direction`` determines the force direction, ``-1`` or for inward and ``+1`` for outward.

 ::

    rhomboid = Rhomboid(pos=[5.0, 5.0, 5.0], a=[1.0, 1.0, 0.0], b=[0.0, 0.0, 1.0], c=[0.0, 1.0, 0.0], direction=1)

creates a rhomboid defined by one corner located at ``[5.0, 5.0, 5.0]`` and three
adjacent edges, defined by the three vectors connecting the corner with its three neighboring corners, ``(1,1,0)`` , ``(0,0,1)`` and ``(0,1,0)``.


:class:`espressomd.shapes.Maze`
    Spherical cavities on a regular grid that are connected by tubes.

The resulting surface is ``nsphere`` spheres of radius ``sphrad`` along each dimension, connected by cylinders of radius ``cylrad``.
The sphere grid have simple cubic symmetry.
The spheres are distributed evenly by dividing the boxl by ``nsphere``.
Dimension of the maze can be controlled by ``dim``: 0 for one dimensional, 1 for two dimensional and 2 for three dimensional maze.


.. figure:: figures/shape-maze.png
   :alt: Example constraint with a ``Maze`` shape.
   :align: center
   :height: 6.00000cm

Pictured is an example constraint with a ``Maze`` shape created with ::

    maze=Maze(cylrad = 2, dim = 2, nsphere = 5, sphrad = 6)
    system.constraints.add(shape=maze, particle_type = 0, penetrable = 1)


:class:`espressomd.shapes.Pore`
    A cylinder with a conical pore between the faces.
  
:todo: `This shape is currently broken. Please do not use.`
    
The pore openings are smoothed with torus segment. The outer radius can be chosen such that it is bigger than the box, to get a wall with a pore. The resulting surface is a cylindrical pore similar to :class:`espressomd.shapes::Cylinder` with a center ``center`` and radius ``radius``.

The ``length`` parameter is half of the cylinder length.
The parameter ``axis`` is a vector along the cylinder axis, which is normalized in the program.
Optionally the outer radius ``outer_rad_left`` and ``outer_rad_right`` of the pore can be specified.
By default these are (numerical) infinity and thus results in an infinite wall with one pore.
The argument radius ``radius`` can be replaced by the argument ``rad_left`` and ``rad_right`` 
to obtain a pore with a conical shape and corresponding opening radii. 
The first radius ``radius_left`` is in the direction opposite to the axis vector.
The same applies for ``outer_radius`` which can be replaced with ``outer_rad_left`` and ``outer_rad_right``.
Per default sharp edges are replaced by circles of unit radius.
The radius of this smoothing can be set with the optional keyword ``smoothing_radius``.

.. figure:: figures/shape-pore1.png
   :alt: Example constraint with a ``Pore`` shape.
   :align: center
   :height: 6.00000cm

Pictured is an example constraint with a ``Pore`` shape created with ::

    pore=Pore(axis = [1,0,0], length = 70, outer_rad_left = 20, outer_rad_right = 30, pos = [50,50,50], rad_left = 10, rad_right = 20, smoothing_radius = 5)
    system.constraints.add(shape=pore, particle_type = 0, penetrable  = 1)

    
:class:`espressomd.shapes.Stomatocyte`
    A stomatocyte.

The resulting surface is a stomatocyte shaped boundary. 
This command should be used with care. 
The position can be any point in the simulation box, 
and the orientation of the (cylindrically symmetric) stomatocyte is given by a vector, 
which points in the direction of the symmetry axis, 
it does not need to be normalized. 
The parameters: ``outer_radius``, ``inner_radius``, and ``layer_width``, specify the shape of the stomatocyte.
Here inappropriate choices of these parameters can yield undersired results. 
The width ``layer_width`` is used as a scaling parameter.
That is, a stomatocyte given by ``outer_radius``:``inner_radius``:``layer_width`` = 7:3:1 
is half the size of the stomatocyte given by 7:3:2. 
Not all choices of the parameters give reasonable values for the shape of the stomatocyte, 
but the combination 7:3:1 is a good point to start from when trying to modify the shape.


.. figure:: figures/shape-stomatocyte1.png
   :alt: Example constraint with a ``Stomatocyte`` shape.
   :align: center
   :height: 6.00000cm

.. figure:: figures/shape-stomatocyte2.png
   :alt: Close-up of the internal ``Stomatocyte`` structure.
   :align: center
   :height: 6.00000cm

   
Pictured is an example constraint with a ``Stomatocyte`` shape (with a closeup of the internal structure) created with ::
  
    stomatocyte=Stomatocyte(inner_radius = 3, outer_radius = 7, orientation_x = 1.0, orientation_y = 0.0,orientation_z = 0.0, position_x = 25, position_y = 25, position_z = 25, layer_width = 3,    direction = 1)
    system.constraints.add(shape=stomatocyte, particle_type = 0, penetrable = 1)

    

:class:`espressomd.shapes.Slitpore`
   Channel-like surface

The resulting surface is T-shape channel that extends in the z-direction.
The cross sectional geometry is depicted in Fig.[fig:slitpore].
It is translationally invariant in y direction.

The region is described as a pore (lower vertical part of the "T"-shape) and a channel (upper horizontal part of the "T"-shape).

.. figure:: figures/slitpore.pdf
   :alt: Schematic for the slitpore shape showing geometrical parameters
   :align: center
   :height: 6.00000cm
   
The parameter ``channel_width`` specifies the distance between the top and the the plateau edge.
The parameter ``pore_length`` specifies the distance between the bottom and the plateau edge.
The parameter ``pore_width`` specifies the distance between the two plateau edges, it is the space between the left and right walls of the pore region.
The parameter ``pore_mouth`` specifies the location (z-coordinate) of the pore opening (centre). It is always centered in the x-direction.

All the edges  are smoothed via the parameters ``upper_smoothing_radius`` (for the concave corner at the edge of the plateau region) and ``lower_smoothing_radius`` (for the convex corner at the bottom of the pore region).
The meaning of the geometrical parameters can be inferred from the shcematic in fig. [fig:slitpore].


.. figure:: figures/shape-slitpore.png
   :alt: Example constraint with a ``Slitpore`` shape.
   :align: center
   :height: 6.00000cm

  
Pictured is an example constraint with a ``Slitpore`` shape created with ::
  
    slitpore=Slitpore(Slitpore(channel_width = 30, lower_smoothing_radius = 3, upper_smoothing_radius = 3, pore_length = 40, pore_mouth = 60, pore_width = 10)
    system.constraints.add(shape=slitpore, particle_type = 0, penetrable = 1)


:class:`espressomd.shapes.SpheroCylinder`
    A capsule, pill, or spherocylinder.
    
The resulting surface is a cylinder capped by hemispheres on both ends.
Similar to `espressomd.shapes::Cylinder`, it is positioned at ``center`` and has a radius ``radius``.
The ``length`` parameter is **half** of the cylinder length, and does not include the contribution from the hemispherical ends.
The ``axis`` parameter is a vector along the cylinder axis, which is normalized in the program.
The direction ``direction`` determines the force direction, ``-1`` or for inward and ``+1`` for outward.


.. figure:: figures/shape-spherocylinder.png
   :alt: Example constraint with a ``SpheroCylinder`` shape.
   :align: center
   :height: 6.00000cm
   
Pictured is an example constraint with a ``SpheroCylinder`` shape created with ::

    spherocylinder = SpheroCylinder(center=[25, 25, 25], axis = [1, 0, 0], direction = 1, radius = 10, length = 30)
    system.constraints.add(shape=spherocylinder, particle_type = 0)


:class:`espressomd.shapes.Hollowcone`
   A hollow cone.

The resulting surface is a section of a hollow cone.
The parameters ``inner_radius`` and ``outer_radius`` specifies the two radii .
The parameter ``opening_angle`` specifies the opening angle of the cone (in radians, between 0 and:math:`\pi/2` ), and thus also determines the length.

The orientation of the (cylindrically symmetric) codne is specified with the parameters ``orientation_x``, ``orientation_y`` and ``orientation_z``. It points in the direction of the symmetry axis, and does not need to be normalized.

The position is specified with ``position_x``, ``position_y`` and ``position_z`` can be any point in the simulation box.

The ``width`` specifies the width.
This shape supports the ``direction`` parameter, +1 the normal points out of the mantel, -1 for when points inward.

.. figure:: figures/shape-hollowcone.png
   :alt:  Example constraint with a  ``Hollowcone`` shape.
   :align: center
   :height: 6.00000cm


Pictured is an example constraint with a ``Hollowcone`` shape created with ::
  
    hollowcone=Hollowcone(HollowCone(inner_radius = 5, outer_radius = 20, opening_angle = np.pi/4.0, orientation_x = 1.0, orientation_y = 0.0, orientation_z = 0.0, position_x = 25, position_y = 25, positi    on_z = 25, width = 2,direction = 1)
    system.constraints.add(shape=hollowcone, particle_type = 0, penetrable = 1)


For the shapes ``wall``; ``sphere``; ``cylinder``; ``rhomboid``; ``maze``; ``pore`` and ``stomacyte``, constraints are able to be penetrated if ``penetrable`` is set to ``True``.
Otherwise, when the ``penetrable`` option is
ignored or is set to `False`, the constraint cannot be violated, i.e. no
particle can go through the constraint surface (|es| will exit if it does).


In variants ``wall``; ``sphere``; ``cylinder``; ``rhomboid`` and ``stomacyte`` it is
also possible to specify a flag indicating if the constraints should be
reflecting. The flags can equal 1 or 2. The flag 1 corresponds to a
reflection process where the normal component of the velocity is
reflected and the tangential component remains unchanged. If the flag is
2, also the tangential component is turned around, so that a bounce back
motion is performed. The second variant is useful for boundaries of DPD.
The reflection property is only activated if an interaction is defined
between a particular particle and the constraint! This will usually be a
lennard-jones interaction with :math:`\epsilon=0`, but finite
interaction range.

..
    .. _Creating a harmonic trap:

    Creating a harmonic trap
    ------------------------

    :todo: `This feature is not yet implemented .`

    Calculates a spring force for all particles, where the equilibrium
    position of the spring is at and its force constant is . A more
    flexible trap can be constructed with constraints, but this one runs on
    the GPU.

.. _Homogeneous Magnetic Field:

Homogeneous Magnetic Field 
--------------------------

:class:`espressomd.Constraints::HomogeneousMagneticField`

This does not define a surface but is based on magnetic dipolar
interaction with an external magnetic field. It applies to all particles
with a dipole moment.

