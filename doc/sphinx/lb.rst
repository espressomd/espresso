.. _Lattice-Boltzmann:

Lattice-Boltzmann
=================

For an implicit treatment of a solvent, |es| allows to couple the molecular
dynamics simulation to a lattice-Boltzmann fluid. The lattice-Boltzmann method (LBM) is a fast, lattice-based method that, in its
"pure" form, allows to calculate fluid flow in different boundary
conditions of arbitrarily complex geometries. Coupled to molecular
dynamics, it allows for the computationally efficient inclusion of
hydrodynamic interactions into the simulation. The focus of the |es| implementation
of the LBM is, of course, the coupling to MD and therefore available
geometries and boundary conditions are somewhat limited in comparison to
"pure" LB codes.

Here we restrict the documentation to the interface. For a more detailed
description of the method, please refer to the literature.

.. note:: Please cite :cite:`arnold13a` (Bibtex key ``arnold13a`` in :file:`doc/sphinx/zrefs.bib`) if you use the LB fluid and :cite:`rohm12a` (Bibtex key ``rohm12a`` in :file:`doc/sphinx/zrefs.bib`) if you use the GPU implementation.

.. _Setting up a LB fluid:

Setting up a LB fluid
---------------------

The following minimal example illustrates how to use the LBM in |es|::

    import espressomd
    system = espressomd.System(box_l=[10, 20, 30])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    lb = espressomd.lb.LBFluid(agrid=1.0, dens=1.0, visc=1.0, tau=0.01)
    system.actors.add(lb)
    system.integrator.run(100)

To use the GPU accelerated variant, replace line 5 in the example above by::

    lb = espressomd.lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, tau=0.01)

.. note:: Feature ``CUDA`` required for GPU accelerated variant

To use the (much faster) GPU implementation of the LBM, use
:class:`espressomd.lb.LBFluidGPU` in place of :class:`espressomd.lb.LBFluid`.
Please note that the GPU implementation uses single precision floating point operations. This decreases the accuracy of calculations compared to the CPU implementation. In particular, due to rounding errors, the fluid density decreases over time, when external forces, coupling to particles, or thermalization is used. The loss of density is on the order of :math:`10^{-12}` per time step.

The command initializes the fluid with a given set of parameters. It is
also possible to change parameters on the fly, but this will only rarely
be done in practice. Before being able to use the LBM, it is necessary
to set up a box of a desired size. The parameter is used to set the
lattice constant of the fluid, so the size of the box in every direction
must be a multiple of ``agrid``.

In the following, we discuss the parameters that can be supplied to the LBM in |es|. The detailed interface definition is available at :class:`espressomd.lb.LBFluid`.

The LB scheme and the MD scheme are not synchronized: In one LB time
step typically several MD steps are performed. This allows to speed up
the simulations and is adjusted with the parameter ``tau``, the LB time step.
The parameters ``dens`` and ``visc`` set up the density and (kinematic) viscosity of the
LB fluid in (usual) MD units. Internally the LB implementation works
with a different set of units: all lengths are expressed in ``agrid``, all times
in ``tau`` and so on.
LB nodes are located at 0.5, 1.5, 2.5, etc.
(in terms of ``agrid``). This has important implications for the location of
hydrodynamic boundaries which are generally considered to be halfway
between two nodes for flat, axis-aligned walls. For more complex boundary geometries, the hydrodynamic boundary location deviates from this midpoint and the deviation decays to first order in ``agrid``.
The LBM should
*not be used as a black box*, but only after a careful check of all
parameters that were applied.

In the following, we describe a number of optional parameters.
Thermalization of the fluid (and particle coupling later on) can be activated by
providing a non-zero value for the parameter ``kT``. Then, a seed has to be provided for
the fluid thermalization::

    lbfluid = espressomd.lb.LBFluid(kT=1.0, seed=134, ...)

The parameter ``ext_force_density`` takes a three dimensional vector as an
array_like of :obj:`float`, representing a homogeneous external body force density in MD
units to be applied to the fluid. The parameter ``bulk_visc`` allows one to
tune the bulk viscosity of the fluid and is given in MD units. In the limit of
low Mach number, the flow does not compress the fluid and the resulting flow
field is therefore independent of the bulk viscosity. It is however known that
the value of the viscosity does affect the quality of the implemented
link-bounce-back method. ``gamma_even`` and ``gamma_odd`` are the relaxation
parameters for the kinetic modes. These fluid parameters do not correspond to
any macroscopic fluid properties, but do influence numerical properties of the
algorithm, such as the magnitude of the error at boundaries. Unless you are an
expert, leave their defaults unchanged. If you do change them, note that they
are to be given in LB units.

Before running a simulation at least the following parameters must be
set up: ``agrid``, ``tau``, ``visc``, ``dens``. For the other parameters, the following are taken: ``bulk_visc=0``, ``gamma_odd=0``, ``gamma_even=0``, ``ext_force_density=[0,0,0]``.

.. _Checkpointing LB:

Checkpointing LB
----------------

::

    lb.save_checkpoint(path, binary)
    lb.load_checkpoint(path, binary)

The first command saves all of the LB fluid nodes' populations to an ascii
(``binary=False``) or binary (``binary=True``) format respectively. The load command
loads the populations from a checkpoint file written with
``lb.save_checkpoint``. In both cases ``path`` specifies the location of the
checkpoint file. This is useful for restarting a simulation either on the same
machine or a different machine. Some care should be taken when using the binary
format as the format of doubles can depend on both the computer being used as
well as the compiler. One thing that one needs to be aware of is that loading
the checkpoint also requires the user to reuse the old forces. This is
necessary since the coupling force between the particles and the fluid has
already been applied to the fluid. Failing to reuse the old forces breaks
momentum conservation, which is in general a problem. It is particularly
problematic for bulk simulations as the system as a whole acquires a drift of
the center of mass, causing errors in the calculation of velocities and
diffusion coefficients. The correct way to restart an LB simulation is to first
load in the particles with the correct forces, and use::

    system.integrator.run(steps=number_of_steps, reuse_forces=True)

upon the first call ``integrator.run``. This causes the
old forces to be reused and thus conserves momentum.

.. _Interpolating velocities:

Interpolating velocities
------------------------

To get interpolated velocity values between lattice nodes, the function::

    lb.get_interpolated_velocity(pos = [1.1,1.2,1.3])
    
with a single position  ``pos`` as an argument can be used. 
For the GPU fluid :class:`espressomd.lb.LBFluidGPU`
also :py:meth:`espressomd.lb.LBFluidGPU.get_interpolated_fluid_velocity_at_positions()`
is available, which expects a numpy array of positions as an argument.

By default, the interpolation is done linearly between the nearest 8 LB nodes,
but for the GPU implementation also a quadratic scheme involving 27 nodes is implemented 
(see eqs. 297 and 301 in :cite:`duenweg08a`). 
You can choose by calling 
one of::

    lb.set_interpolation_order('linear')
    lb.set_interpolation_order('quadratic')
    
A note on boundaries:
both interpolation schemes don't take into account the physical location of the boundaries
(e.g. in the middle between two nodes for a planar wall) but will use the boundary node slip velocity 
at the node position. This means that every interpolation involving at least one
boundary node will introduce an error.

.. _Coupling LB to a MD simulation:

Coupling LB to a MD simulation
------------------------------

MD particles can be coupled to a LB fluid through frictional coupling. The friction force

  .. math:: F_{i,\text{frict}} = - \gamma (v_i(t)-u(x_i(t),t))

depends on the particle velocity :math:`v` and the fluid velocity :math:`u`. It acts both
on the particle and the fluid (in opposite direction). Because the fluid is also affected,
multiple particles can interact via hydrodynamic interactions. As friction in molecular systems is
accompanied by fluctuations, the particle-fluid coupling has to be activated through
the :ref:`LB thermostat` (See more detailed description there). A short example is::

    system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

where ``lbf`` is an instance of either :class:`espressomd.lb.LBFluid` or :class:`espressomd.lb.LBFluidGPU`, 
``gamma`` the friction coefficient and ``seed`` the seed for the random number generator involved 
in the thermalization.


.. _Reading and setting properties of single lattice nodes:

Reading and setting properties of single lattice nodes
------------------------------------------------------

Appending three indices to the ``lb`` object returns an object that represents the selected LB grid node and allows one to access all of its properties::

    lb[x, y, z].density     # fluid density (one scalar for LB and CUDA)
    lb[x, y, z].velocity    # fluid velocity (a numpy array of three floats)
    lb[x, y, z].stress      # fluid pressure tensor (a symmetric 3x3 numpy array of floats)
    lb[x, y, z].stress_neq  # nonequilbrium part of the pressure tensor (as above)
    lb[x, y, z].boundary    # flag indicating whether the node is fluid or boundary (fluid: boundary=0, boundary: boundary != 0)
    lb[x, y, z].population  # 19 LB populations (a numpy array of 19 floats, check order from the source code)

All of these properties can be read and used in further calculations. Only the property ``population`` can be modified. The indices ``x,y,z`` are integers and enumerate the LB nodes in the three directions, starts with 0. To modify ``boundary``, refer to :ref:`Setting up boundary conditions`.

Examples::

    print(lb[0, 0, 0].velocity)

    lb[0, 0, 0].density = 1.2

The first line prints the fluid velocity at node 0 0 0 to the screen. The second line sets this fluid node's density to the value ``1.2``.

.. _Removing total fluid momentum:

Removing total fluid momentum
-----------------------------

.. note:: Only available for ``CUDA``

Some simulations require the net momentum of the system to vanish. Even if the
physics of the system fulfills this condition, numerical errors can introduce
drift. To remove the momentum in the fluid call::

    lb.remove_momentum()

.. _Output for visualization:

Output for visualization
------------------------

|es| implements a number of commands to output fluid field data of the whole fluid into a file at once. ::

    lb.print_vtk_velocity(path)
    lb.print_vtk_boundary(path)
    lb.print_velocity(path)
    lb.print_boundary(path)

Currently supported fluid properties are the velocity, and boundary flag in ASCII VTK as well as Gnuplot compatible ASCII output.

The VTK format is readable by visualization software such as ParaView [1]_
or Mayavi2 [2]_. If you plan to use ParaView for visualization, note that also the particle
positions can be exported using the VTK format (see :meth:`~espressomd.particle_data.ParticleList.writevtk`).

The variant

::

   lb.print_vtk_velocity(path, bb1, bb2)

allows you to only output part of the flow field by specifying an axis aligned
bounding box through the coordinates ``bb1`` and ``bb1`` (lists of three ints) of two of its corners. This
bounding box can be used to output a slice of the flow field. As an
example, executing

::

    lb.print_vtk_velocity(path, [0, 0, 5], [10, 10, 5])

will output the cross-section of the velocity field in a plane
perpendicular to the :math:`z`-axis at :math:`z = 5` (assuming the box
size is 10 in the :math:`x`- and :math:`y`-direction).

.. If the bicomponent fluid is used, two filenames have to be supplied when exporting the density field, to save both components.


.. _Choosing between the GPU and CPU implementations:

Choosing between the GPU and CPU implementations
------------------------------------------------

.. note:: Feature ``CUDA`` required

|es| contains an implementation of the LBM for NVIDIA
GPUs using the CUDA framework. On CUDA-supporting machines this can be
activated by compiling with the feature ``CUDA``. Within the
Python script, the :class:`~espressomd.lb.LBFluid` object can be substituted with the :class:`~espressomd.lb.LBFluidGPU` object to switch from CPU based to GPU based execution. For further
information on CUDA support see section :ref:`GPU Acceleration with CUDA`.

The following minimal example demonstrates how to use the GPU implementation of the LBM in analogy to the example for the CPU given in section :ref:`Setting up a LB fluid`::

    import espressomd
    system = espressomd.System(box_l=[10, 20, 30])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    lb = espressomd.lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, tau=0.01)
    system.actors.add(lb)
    system.integrator.run(100)

For boundary conditions analogous to the CPU
implementation, the feature ``LB_BOUNDARIES_GPU`` has to be activated.
The feature ``CUDA`` allows the use of Lees-Edwards boundary conditions. Our implementation follows the paper of :cite:`wagner02`. Note, that there is no extra python interface for the use of Lees-Edwards boundary conditions with the LB algorithm. All information are rather internally derived from the set of the Lees-Edwards offset in the system class. For further information Lees-Edwards boundary conditions please refer to section :ref:`Lees-Edwards boundary conditions`

.. _Electrohydrodynamics:

Electrohydrodynamics
--------------------

        .. note::
           This needs the feature ``LB_ELECTROHYDRODYNAMICS``.

If the feature is activated, the lattice-Boltzmann code can be
used to implicitly model surrounding salt ions in an external electric
field by having the charged particles create flow.

For that to work, you need to set the electrophoretic mobility
(multiplied by the external :math:`E`-field) :math:`\mu E` on the
particles that should be subject to the field. This effectively acts
as a velocity offset between the particle and the LB fluid.

For more information on this method and how it works, read the
publication :cite:`hickey10a`.


.. _Using shapes as lattice-Boltzmann boundary:

Using shapes as lattice-Boltzmann boundary
------------------------------------------

.. note::
    Feature ``LB_BOUNDARIES`` required

Lattice-Boltzmann boundaries are implemented in the module
:mod:`espressomd.lbboundaries`. You might want to take a look
at the classes :class:`espressomd.lbboundaries.LBBoundary`
and :class:`espressomd.lbboundaries.LBBoundaries` for more information.

Adding a shape-based boundary is straightforward::

    lbb = espressomd.lbboundaries.LBBoundary(shape=my_shape, velocity=[0, 0, 0])
    system.lbboundaries.add(lbb)

or::

    lbb = espressomd.lbboundaries.LBBoundary()
    lbb.shape = my_shape
    lbb.velocity = [0, 0, 0]
    system.lbboundaries.add(lbb)

.. _Minimal usage example:

Minimal usage example
~~~~~~~~~~~~~~~~~~~~~

.. note:: Feature ``LB_BOUNDARIES`` or ``LB_BOUNDARIES_GPU`` required

In order to add a wall as boundary for a lattice-Boltzmann fluid
you could do the following::

    wall = espressomd.shapes.Wall(dist=5, normal=[1, 0, 0])
    lbb = espressomd.lbboundaries.LBBoundary(shape=wall, velocity=[0, 0, 0])
    system.lbboundaries.add(lbb)

.. _Setting up boundary conditions:

Setting up boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example sets up a system consisting of a spherical boundary in the center of the simulation box acting as a no-slip boundary for the LB fluid that is driven by 4 walls with a slip velocity::

    from espressomd import System, lb, lbboundaries, shapes

    system = System(box_l=[64, 64, 64])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    lb = lb.LBFluid(agrid=1.0, dens=1.0, visc=1.0, tau=0.01)
    system.actors.add(lb)

    v = [0, 0, 0.01]  # the boundary slip
    walls = [None] * 4

    wall_shape = shapes.Wall(normal=[1, 0, 0], dist=1)
    walls[0] = lbboundaries.LBBoundary(shape=wall_shape, velocity=v)

    wall_shape = shapes.Wall(normal=[-1, 0, 0], dist=-63)
    walls[1] = lbboundaries.LBBoundary(shape=wall_shape, velocity=v)

    wall_shape = shapes.Wall(normal=[0, 1, 0], dist=1)
    walls[2] = lbboundaries.LBBoundary(shape=wall_shape, velocity=v)

    wall_shape = shapes.Wall(normal=[0, -1, 0], dist=-63)
    walls[3] = lbboundaries.LBBoundary(shape=wall_shape, velocity=v)

    for wall in walls:
        system.lbboundaries.add(wall)

    sphere_shape = shapes.Sphere(radius=5.5, center=[33, 33, 33], direction=1)
    sphere = lbboundaries.LBBoundary(shape=sphere_shape)
    system.lbboundaries.add(sphere)

    system.integrator.run(4000)

    print(sphere.get_force())

After integrating the system for a sufficient time to reach the steady state, the hydrodynamic drag force exerted on the sphere is evaluated.

The LB boundaries use the same :mod:`~espressomd.shapes` objects to specify their geometry as :mod:`~espressomd.constraints` do for particles. This allows the user to quickly set up a system with boundary conditions that simultaneously act on the fluid and particles. For a complete description of all available shapes, refer to :mod:`espressomd.shapes`.

Intersecting boundaries are in principle possible but must be treated
with care. In the current implementation, all nodes that are
within at least one boundary are treated as boundary nodes.

Currently, only the so-called "link-bounce-back" algorithm for wall
nodes is available. This creates a boundary that is located
approximately midway between the lattice nodes, so in the above example ``wall[0]``
corresponds to a boundary at :math:`x=1.5`. Note that the
location of the boundary is unfortunately not entirely independent of
the viscosity. This can be seen when using the sample script with a high
viscosity.

The bounce back boundary conditions permit it to set the velocity at the boundary
to a nonzero value via the ``v`` property of an ``LBBoundary`` object. This allows to create shear flow and boundaries
moving relative to each other. The velocity boundary conditions are
implemented according to :cite:`succi01a` eq. 12.58. Using
this implementation as a blueprint for the boundary treatment, an
implementation of the Ladd-Coupling should be relatively
straightforward. The ``LBBoundary`` object furthermore possesses a property ``force``, which keeps track of the hydrodynamic drag force exerted onto the boundary by the moving fluid.


.. [1]
   http://www.paraview.org/

.. [2]
   http://code.enthought.com/projects/mayavi/
