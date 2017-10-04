.. _Setting up shapes:

=================
Setting up shapes
=================

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

Creating different shapes
=========================
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


Using shapes as constraints
===========================

.. note::
    `Feature CONSTRAINTS required`

Adding shape based constraints to the system
--------------------------------------------


Usually you want to use constraints based on a shape.
The module :mod:`espressomd.constraints` provides the class
:class:`espressomd.constraints.ShapeBasedConstraint`::

    shape_constraint = espressomd.constraints.ShapeBasedConstraint(shape=my_shape)

In order to add the constraint to to the system
invoke the :meth:`espressomd.constraints.add` method::

    system.constraints.add(shape_constraint)


Usage example
-------------

If we wanted to add a non-penetrable pore constraint to our simulation,
we could do the following::

    pore = espressomd.shapes.Pore(axis=[1,0,0], length=2, pos=[15,15,15], smoothing_radius=0.5)
    pore_constraint = espressomd.constraints.ShapeBasedConstraint(shape=pore, penetrable=0, particle_type=1)
    system.constraints.add(pore_constraint)

Interactions between the pore and other particles are then defined
as usual (:ref:`Setting up interactions`).


Using shapes as Lattice-Boltzmann Boundary
==========================================

.. note::
    `Feature LB_BOUNDARIES required`

Adding shape based Lattice-Boltzmann Boundary
---------------------------------------------

Lattice-Boltzmann boundaries are implemented in the module
:mod:`espressomd.lbboundaries`. You might want to take a look
at the classes :class:`espressomd.lbboundaries.LBBoundary`
and :class:`espressomd.lbboundaries.LBBoundaries` for more information.

Adding a shape based boundary is straightforward::

    lbb = espressomd.lbboundaries.LBBoundary(shape=my_shape, velocity=[0,0,0])
    system.lbboundaries.add(lbb)

or::

    lbb = espressomd.lbboundaries.LBBoundary()
    lbb.shape = my_shape
    lbb.velocity = [0,0,0]
    system.lbboundaries.add(lbb)

Usage example
-------------

In order to add a wall as boundary for a Lattice-Boltzmann fluid
you could do the following::

    wall = espressomd.shapes.Wall(dist=5, normal=[1,0,0])
    lbb = espressomd.lbboundaries.LBBoundary(shape=wall, velocity=[0,0,0])
    system.lbboundaries.add(lbb)
