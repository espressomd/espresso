.. _Lattice-Boltzmann:

Lattice-Boltzmann
=================

For an implicit treatment of a solvent, |es| can couple the molecular
dynamics simulation to a lattice-Boltzmann fluid. The lattice-Boltzmann
method (LBM) is a fast, lattice-based method that, in its "pure" form,
allows to calculate fluid flow in different boundary conditions of
arbitrarily complex geometries. Coupled to molecular dynamics,
it allows for the computationally efficient inclusion of hydrodynamic
interactions into the simulation. The focus of the |es| implementation
of the LBM is, of course, the coupling to MD and therefore available
geometries and boundary conditions are somewhat limited in comparison to
"pure" LB codes.

Here we restrict the documentation to the interface. For a more detailed
description of the method, please refer to the literature.

.. note::
    Please cite :cite:t:`godenschwager13a` and :cite:t:`bauer21a` (BibTeX keys
    ``godenschwager13a`` and ``bauer21a`` in :file:`doc/bibliography.bib`) if
    you use the LB fluid. When generating your own kernels with pystencils and
    lbmpy, please also cite :cite:t:`bauer19a` and :cite:t:`bauer21b` (BibTeX
    key ``bauer19a`` resp. ``bauer21b`` in :file:`doc/bibliography.bib`).

.. note::

    Requires external feature ``WALBERLA``, enabled with the CMake option
    ``-D ESPRESSO_BUILD_WITH_WALBERLA=ON``.

.. _Setting up a LB fluid:

Setting up a LB fluid
---------------------

The following minimal example illustrates how to use the LBM in |es|::

    import espressomd
    import espressomd.lb
    system = espressomd.System(box_l=[10, 20, 30])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    lb = espressomd.lb.LBFluidWalberla(agrid=1.0, density=1.0, viscosity=1.0, tau=0.01)
    system.actors.add(lb)
    system.integrator.run(100)

To use the GPU-accelerated variant, replace line 6 in the example above by::

    lb = espressomd.lb.LBFluidWalberlaGPU(agrid=1.0, density=1.0, viscosity=1.0, tau=0.01)

.. note:: Feature ``CUDA`` required for the GPU-accelerated variant

To use the (much faster) GPU implementation of the LBM, use
:class:`~espressomd.lb.LBFluidWalberlaGPU` in place of :class:`~espressomd.lb.LBFluidWalberla`.
Please note that the GPU implementation uses single precision floating point operations.
This decreases the accuracy of calculations compared to the CPU implementation.
In particular, due to rounding errors, the fluid density decreases over time,
when external forces, coupling to particles, or thermalization is used.
The loss of density is on the order of :math:`10^{-12}` per time step.

The command initializes the fluid with a given set of parameters. It is
also possible to change parameters on the fly, but this will only rarely
be done in practice. Before being able to use the LBM, it is necessary
to set up a box of a desired size. The parameter is used to set the
lattice constant of the fluid, so the size of the box in every direction
must be a multiple of ``agrid``.

In the following, we discuss the parameters that can be supplied to the LBM in |es|.
The detailed interface definition is available at :class:`~espressomd.lb.LBFluidWalberla`.

The LB scheme and the MD scheme are not synchronized: In one LB time
step typically several MD steps are performed. This allows to speed up
the simulations and is adjusted with the parameter ``tau``, the LB time step.
The parameters ``density`` and ``viscosity`` set up the density and (kinematic) viscosity of the
LB fluid in (usual) MD units. Internally the LB implementation works
with a different set of units: all lengths are expressed in ``agrid``, all times
in ``tau`` and so on.
LB nodes are located at 0.5, 1.5, 2.5, etc.
(in terms of ``agrid``). This has important implications for the location of
hydrodynamic boundaries which are generally considered to be halfway
between two nodes for flat, axis-aligned walls. For more complex boundary geometries,
the hydrodynamic boundary location deviates from this midpoint and the deviation
decays to first order in ``agrid``. The LBM should
*not be used as a black box*, but only after a careful check of all
parameters that were applied.

In the following, we describe a number of optional parameters.
Thermalization of the fluid (and particle coupling later on) can be activated by
providing a non-zero value for the parameter ``kT``. Then, a seed has to be provided for
the fluid thermalization::

    lb = espressomd.lb.LBFluidWalberla(kT=1.0, seed=134, ...)

The parameter ``ext_force_density`` takes a three dimensional vector as an
array_like of :obj:`float`, representing a homogeneous external body force density in MD
units to be applied to the fluid.

Before running a simulation at least the following parameters must be
set up: ``agrid``, ``tau``, ``viscosity``, ``density``.

Performance considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The CPU implementation of the LB has an extra flag ``single_precision`` to
use single-precision floating point values. These are approximately 10%
faster than double-precision, at the cost of a small loss in precision.

To enable vectorization, run ``cmake . -DWALBERLA_USE_AVX=ON``.
An AVX2-capable microprocessor is required. Currently only works
for double-precision kernels.

.. _Checkpointing LB:

Checkpointing
-------------

::

    lb.save_checkpoint(path, binary)
    lb.load_checkpoint(path, binary)

The first command saves all of the LB fluid nodes' populations to an ASCII
(``binary=False``) or binary (``binary=True``) format respectively.
The second command loads the LB fluid nodes' populations.
In both cases ``path`` specifies the location of the
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

    lb.get_interpolated_velocity(pos=[1.1, 1.2, 1.3])

with a single position  ``pos`` as an argument can be used.

The interpolation is done linearly between the nearest 8 LB nodes.

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
the :ref:`LB thermostat` (see more detailed description there). A short example is::

    system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

where ``lbf`` is an instance of either :class:`~espressomd.lb.LBFluidWalberla` or
:class:`~espressomd.lb.LBFluidWalberlaGPU`, ``gamma`` the friction coefficient and
``seed`` the seed for the random number generator involved
in the thermalization.

.. _LB and LEbc:

LB and LEbc
^^^^^^^^^^^

:ref:`Lees-Edwards boundary conditions` (LEbc) are supported by both
LB implementations, which follow the derivation in :cite:`wagner02a`.
Note, that there is no extra python interface for the use of LEbc
with the LB algorithm: all the necessary information is internally
derived from the currently active MD LEbc protocol in
``system.lees_edwards.protocol``.
Therefore, the MD LEbc must be set before the LB actor is instantiated.
Use the :class:`~espressomd.lees_edwards.Off` if the system should have
no shearing initially; this action will initialize the shear axes, and
when the LB actor is instantiated, the Lees-Edwards collision kernels
will be used instead of the default ones.

.. note::

    At the moment, LB only supports the case ``shear_plane_normal="y"``.

.. _Reading and setting properties of single lattice nodes:

Reading and setting properties of single lattice nodes
------------------------------------------------------

Appending three indices to the ``lb`` object returns an object that represents
the selected LB grid node and allows one to access all of its properties::

    lb[x, y, z].density              # fluid density (one scalar for LB and CUDA)
    lb[x, y, z].velocity             # fluid velocity (a numpy array of three floats)
    lb[x, y, z].pressure_tensor      # fluid pressure tensor (a symmetric 3x3 numpy array of floats)
    lb[x, y, z].is_boundary          # flag indicating whether the node is fluid or boundary (fluid: boundary=0, boundary: boundary != 1)
    lb[x, y, z].population           # 19 LB populations (a numpy array of 19 floats, check order from the source code)

All of these properties can be read and used in further calculations.
Only the property ``population`` can be modified. The indices ``x, y, z``
are integers and enumerate the LB nodes in the three Cartesian directions,
starting at 0. To modify ``is_boundary``, refer to :ref:`Setting up LB boundary conditions`.

Example::

    print(lb[0, 0, 0].velocity)
    lb[0, 0, 0].density = 1.2

The first line prints the fluid velocity at node (0 0 0) to the screen.
The second line sets this fluid node's density to the value ``1.2``.
Use negative indices to get nodes starting from the end of the lattice.

The nodes can be read and modified using slices. Example::

    print(lb[0:4:2, 0:2, 0].velocity)
    lb[0:4:2, 0:2, 0].density = [[[1.1], [1.2]], [[1.3], [1.4]]]

The first line prints an array of shape (2, 2, 1, 3) with the velocities
of nodes (0 0 0), (0 1 0), (2 0 0), (2 1 0). The second line updates
these nodes with densities ranging from 1.1 to 1.4. You can set either
a value that matches the length of the slice (which sets each node
individually), or a single value that will be copied to every node
(e.g. a scalar for density, or an array of length 3 for the velocity).

.. _LB VTK output:

VTK output
----------

The waLBerla library implements a globally-accessible VTK registry.
A VTK stream can be attached to a LB actor to periodically write
one or multiple fluid field data into a single file using
:class:`~espressomd.lb.VTKOutput`::

    vtk_obs = ["density", "velocity_vector"]
    # create a VTK callback that automatically writes every 10 LB steps
    lb_vtk = espressomd.lb.VTKOutput(
        lb_fluid=lb, identifier="lb_vtk_automatic", observables=vtk_obs,
        delta_N=10)
    self.system.integrator.run(100)
    # can be deactivated
    lb_vtk.disable()
    self.system.integrator.run(10)
    lb_vtk.enable()
    # create a VTK callback that writes only when explicitly called
    lb_vtk_on_demand = espressomd.lb.VTKOutput(
        lb_fluid=lb, identifier="lb_vtk_now", observables=vtk_obs)
    lb_vtk_on_demand.write()

Currently supported fluid properties are the density, velocity vector
and pressure tensor. By default, the properties of the current state
of the fluid are written to disk on demand. To add a stream that writes
to disk continuously, use the optional argument ``delta_N`` to indicate
the level of subsampling. Such a stream can be deactivated.

The VTK format is readable by visualization software such as ParaView [1]_
or Mayavi2 [2]_. If you plan to use ParaView for visualization, note that also the particle
positions can be exported using the VTK format (see :meth:`~espressomd.particle_data.ParticleList.writevtk`).

.. _Choosing between the GPU and CPU implementations:

Choosing between the GPU and CPU implementations
------------------------------------------------

|es| contains an implementation of the LBM for NVIDIA
GPUs using the CUDA framework. On CUDA-supporting machines this can be
activated by compiling with the feature ``CUDA``. Within the
Python script, the :class:`~espressomd.lb.LBFluidWalberla` object can be substituted
with the :class:`~espressomd.lb.LBFluidWalberlaGPU` object to switch from CPU based
to GPU based execution. For further
information on CUDA support see section :ref:`CUDA acceleration`.

The following minimal example demonstrates how to use the GPU implementation
of the LBM in analogy to the example for the CPU given in section
:ref:`Setting up a LB fluid`::

    import espressomd
    system = espressomd.System(box_l=[10, 20, 30])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    lb = espressomd.lb.LBFluidWalberlaGPU(agrid=1.0, density=1.0, viscosity=1.0, tau=0.01)
    system.actors.add(lb)
    system.integrator.run(100)

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
publication :cite:t:`hickey10a`.

.. _Setting up LB boundary conditions:

Setting up boundary conditions
------------------------------

Currently, only the so-called "link-bounce-back" algorithm for boundary
nodes is available. This creates a boundary that is located
approximately midway between lattice nodes. With no-slip boundary conditions,
populations are reflected back. With slip velocities, the reflection is
followed by a velocity interpolation. This allows to create shear flow and
boundaries "moving" relative to each other.

Under the hood, a boundary field is added to the blockforest, which contains
pre-calculated information for the reflection and interpolation operations.

.. _Per-node LB boundary conditions:

Per-node boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can set (or update) the slip velocity of individual nodes::

    import espressomd.lb
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    lbf = espressomd.lb.LBFluidWalberla(agrid=0.5, density=1.0, viscosity=1.0, tau=0.01)
    system.actors.add(lbf)
    # make one node a boundary node with a slip velocity
    lbf[0, 0, 0].boundary = espressomd.lb.VelocityBounceBack([0, 0, 1])
    # update node for no-slip boundary conditions
    lbf[0, 0, 0].boundary = espressomd.lb.VelocityBounceBack([0, 0, 0])
    # remove boundary conditions
    lbf[0, 0, 0].boundary = None

.. _Shape-based LB boundary conditions:

Shape-based boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding a shape-based boundary is straightforward::

    import espressomd.lb
    import espressomd.shapes
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    lbf = espressomd.lb.LBFluidWalberla(agrid=0.5, density=1.0, viscosity=1.0, tau=0.01)
    system.actors.add(lbf)
    # set up shear flow between two sliding walls
    wall1 = espressomd.shapes.Wall(normal=[+1., 0., 0.], dist=2.5)
    lbf.add_boundary_from_shape(shape=wall1, velocity=[0., +0.05, 0.])
    wall2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-(system.box_l[0] - 2.5))
    lbf.add_boundary_from_shape(shape=wall2, velocity=[0., -0.05, 0.])

The ``velocity`` argument is optional, in which case the no-slip boundary
conditions are used. For a position-dependent slip velocity, the argument
to ``velocity`` must be a 4D grid (the first three dimensions must match
the LB grid shape, the fourth dimension has size 3 for the velocity).

The LB boundaries use the same :mod:`~espressomd.shapes` objects to specify
their geometry as :mod:`~espressomd.constraints` do for particles.
This allows the user to quickly set up a system with boundary conditions
that simultaneously act on the fluid and particles. For a complete
description of all available shapes, refer to :mod:`espressomd.shapes`.


.. [1]
   https://www.paraview.org/

.. [2]
   http://code.enthought.com/projects/mayavi/
