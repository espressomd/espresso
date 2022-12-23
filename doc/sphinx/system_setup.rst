.. _Setting up the system:

Setting up the system
=====================

.. _Setting global variables:

Setting global variables
------------------------

The global system variables are controlled via the Python
:class:`espressomd.system.System` class::

    import espressomd
    system = espressomd.System(box_l=[10., 10., 10.])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

This code creates a system with a cubic unit cell of length 10 in simulation
units, and sets the skin to 0.4 simulation units and the time step to 0.01
simulation units.

Some variables belong to the ``System`` class and will be explained in the list
below, while other variables such as the skin belong to objects that are
attached to the class, and will be explain in subsequent sections and chapters.
Note that for many vectorial properties, e.g. ``box_l`` and ``periodicity``,
component-wise manipulation
like ``system.box_l[0] = 1`` or in-place operators like ``+=`` or ``*=`` are not
allowed and result in an error. This behavior is inherited, so the same applies
to ``a`` after ``a = system.box_l``. If you want to use a vectorial property
for further calculations, you should explicitly make a copy e.g. via
``a = numpy.copy(system.box_l)``.

* :py:attr:`~espressomd.system.System.box_l`

  Simulation box lengths of the cuboid box used by |es|.
  Note that if you change the box length during the simulation, the folded
  particle coordinates will remain the same, i.e., the particle stay in
  the same image box, but at the same relative position in their image
  box. If you want to scale the positions, use the command
  :py:meth:`~espressomd.system.System.change_volume_and_rescale_particles`.

* :py:attr:`~espressomd.system.System.periodicity`

  Specifies periodicity for the three directions. |es| can be instructed
  to treat some dimensions as non-periodic. By default |es| assumes periodicity in
  all directions which equals setting this variable to ``[True, True, True]``.
  A dimension is specified as non-periodic via setting the periodicity
  variable for this dimension to ``False``. E.g. Periodicity only in z-direction
  is obtained by ``[False, False, True]``. Caveat: Be aware of the fact that making a
  dimension non-periodic does not hinder particles from leaving the box in
  this direction; in this case, shape-based constraints can be used to keep
  particles in the simulation box. For more details, see :ref:`Boundary conditions`.

* :py:attr:`~espressomd.system.System.time_step`

  Time step for MD integration.

* :py:attr:`~espressomd.system.System.time`

  The simulation time.

* :py:attr:`~espressomd.system.System.min_global_cut`

  Minimal total cutoff for real space. Effectively, this plus the
  :py:attr:`~espressomd.cell_system.CellSystem.skin` is the minimally possible
  cell size. |es| typically determines this value automatically, but some
  algorithms, virtual sites, require you to specify it manually.

* :py:attr:`~espressomd.system.System.max_cut_bonded`

  *read-only* Maximal cutoff of bonded interactions.

* :py:attr:`~espressomd.system.System.max_cut_nonbonded`

  *read-only* Maximal cutoff of non-bonded interactions.

.. _Accessing module states:

Accessing module states
~~~~~~~~~~~~~~~~~~~~~~~

Some variables like or are no longer directly available as attributes.
In these cases they can be easily derived from the corresponding Python
objects like::

    n_part = len(system.part)

or by calling the corresponding ``get_state()`` methods like::

    temperature = system.thermostat.get_state()[0]['kT']
    gamma = system.thermostat.get_state()[0]['gamma']
    gamma_rot = system.thermostat.get_state()[0]['gamma_rotation']

.. _Simulation box:

Simulation box
--------------

.. _Boundary conditions:

Boundary conditions
~~~~~~~~~~~~~~~~~~~

.. _Periodic boundaries:

Periodic boundaries
^^^^^^^^^^^^^^^^^^^

With periodic boundary conditions, particles interact with periodic
images of all particles in the system. This is the default behavior.
When particles cross a box boundary, their position are folded and
their image box counter are incremented.

From the Python interface, the folded position is accessed with
:attr:`~espressomd.particle_data.ParticleHandle.pos_folded` and the image
box counter with :attr:`~espressomd.particle_data.ParticleHandle.image_box`.
Note that :attr:`~espressomd.particle_data.ParticleHandle.pos` gives the
unfolded particle position.

Example::

    import espressomd
    system = espressomd.System(box_l=[5.0, 5.0, 5.0], periodicity=[True, True, True])
    system.time_step = 0.1
    system.cell_system.skin = 0.0
    p = system.part.add(pos=[4.9, 0.0, 0.0], v=[0.1, 0.0, 0.0])
    system.integrator.run(20)
    print(f"pos        = {p.pos}")
    print(f"pos_folded = {p.pos_folded}")
    print(f"image_box  = {p.image_box}")

Output:

.. code-block:: none

    pos        = [5.1 0.  0. ]
    pos_folded = [0.1 0.  0. ]
    image_box  = [1 0 0]

.. _Open boundaries:

Open boundaries
^^^^^^^^^^^^^^^

With open boundaries, particles can leave the simulation box.
What happens in this case depends on which algorithm is used.
Some algorithms may require open boundaries,
such as :ref:`Stokesian Dynamics`.

Example::

    import espressomd
    system = espressomd.System(box_l=[5.0, 5.0, 5.0], periodicity=[False, False, False])
    system.time_step = 0.1
    system.cell_system.skin = 0.0
    p = system.part.add(pos=[4.9, 0.0, 0.0], v=[0.1, 0.0, 0.0])
    system.integrator.run(20)
    print(f"pos        = {p.pos}")
    print(f"pos_folded = {p.pos_folded}")
    print(f"image_box  = {p.image_box}")

Output:

.. code-block:: none

    pos        = [5.1 0.  0. ]
    pos_folded = [5.1 0.  0. ]
    image_box  = [0 0 0]

.. _Lees-Edwards boundary conditions:

Lees--Edwards boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lees--Edwards boundary conditions (LEbc) are special periodic boundary
conditions to simulate systems under shear stress :cite:`lees72a`.
Periodic images of particles across the shear boundary appear with a
time-dependent position offset. When a particle crosses the shear boundary,
it appears to the opposite side of the simulation box with a position offset
and a shear velocity :cite:`bindgen21a`.

LEbc require a fully periodic system and are configured with
:class:`~espressomd.lees_edwards.LinearShear` and
:class:`~espressomd.lees_edwards.OscillatoryShear`.
To temporarily disable LEbc, use :class:`~espressomd.lees_edwards.Off`.
To completely disable LEbc and reinitialize the box geometry, do
``system.lees_edwards.protocol = None``.

Example::

    import espressomd
    import espressomd.lees_edwards
    system = espressomd.System(box_l=[5.0, 5.0, 5.0])
    system.time_step = 0.1
    system.cell_system.skin = 0.0
    system.cell_system.set_n_square(use_verlet_lists=True)
    le_protocol = espressomd.lees_edwards.LinearShear(
        shear_velocity=-0.1, initial_pos_offset=0.0, time_0=-0.1)
    system.lees_edwards.set_boundary_conditions(
        shear_direction="y", # shear along y-axis
        shear_plane_normal="x", # shift when crossing the x-boundary
        protocol=le_protocol)
    p = system.part.add(pos=[4.9, 0.0, 0.0], v=[0.1, 0.0, 0.0])
    system.integrator.run(20)
    print(f"pos        = {p.pos}")
    print(f"pos_folded = {p.pos_folded}")
    print(f"image_box  = {p.image_box}")
    print(f"velocity   = {p.v}")

Output:

.. code-block:: none

    pos        = [5.1 0.2 0. ]
    pos_folded = [0.1 0.2 0. ]
    image_box  = [1 0 0]
    velocity   = [0.1 0.1 0. ]

Particles inserted outside the box boundaries will be wrapped around
using the normal periodic boundary rules, i.e. they will not be sheared,
even though their :attr:`~espressomd.particle_data.ParticleHandle.image_box`
is *not* zero.

Once a valid tuple ``(shear_direction, shear_plane_normal, protocol)`` has been
set via :meth:`~espressomd.lees_edwards.LeesEdwards.set_boundary_conditions`,
one can update the protocol via a simple assignment of the form
``system.lees_edwards.protocol = new_le_protocol``, in which case
the shear direction and shear normal are left unchanged. The method
:meth:`~espressomd.lees_edwards.LeesEdwards.set_boundary_conditions`
is the only way to modify the shear direction and shear normal.


.. _Cell systems:

Cell systems
~~~~~~~~~~~~

This section deals with the flexible particle data organization of |es|.
|es| is able to change the organization of the particles in the computer
memory to accommodate for the needs of the algorithms being used.
For details on the internal organization,
refer to section :ref:`Internal particle organization`.

.. _Global properties:

Global properties
^^^^^^^^^^^^^^^^^

The properties of the cell system can be accessed via the system
:class:`~espressomd.system.System.cell_system` attribute:

* :py:attr:`~espressomd.cell_system.CellSystem.node_grid`

  3D node grid for real space domain decomposition (optional, if
  unset an optimal partition is chosen automatically). The domain decomposition
  can be visualized with :file:`samples/visualization_cellsystem.py`.

* :py:attr:`~espressomd.cell_system.CellSystem.skin`

  Skin for the Verlet list. This value has to be set, otherwise the simulation will not start.

Details about the cell system can be obtained by
:meth:`get_state() <espressomd.cell_system.CellSystem.get_state>`:

* ``cell_grid``       Dimension of the inner cell grid (only for regular decomposition).
* ``cell_size``       Box-length of a cell (only for regular decomposition).
* ``n_nodes``         Number of MPI nodes.
* ``node_grid``       MPI domain partition.
* ``type``            The current type of the cell system.
* ``skin``            Verlet list skin.
* ``verlet_reuse``    Average number of integration steps the Verlet list is re-used.

.. _Regular decomposition:

Regular decomposition
^^^^^^^^^^^^^^^^^^^^^

Invoking :py:meth:`~espressomd.cell_system.CellSystem.set_regular_decomposition`
selects the regular decomposition cell scheme, using Verlet lists for the
calculation of the interactions. If you specify ``use_verlet_lists=False``,
only the regular decomposition is used, but not the Verlet lists. ::

    import espressomd
    system = espressomd.System(box_l=[1, 1, 1])
    system.cell_system.set_regular_decomposition(use_verlet_lists=True)

The regular decomposition cellsystem is the default system and suits most
applications with short ranged interactions. The particles are divided
up spatially into small compartments, the cells, such that the cell size
is larger than the maximal interaction range. In this case interactions
only occur between particles in adjacent cells. Since the interaction
range should be much smaller than the total system size, leaving out all
interactions between non-adjacent cells can mean a tremendous speed-up.
Moreover, since for constant interaction range, the number of particles
in a cell depends only on the density. The number of interactions is
therefore of the order :math:`N` instead of order :math:`N^2` if one has to
calculate all pair interactions.

With this scheme, there must be at least two cells per direction,
and at most 32 cells per direction for a cubic box geometry.
The number of cells per direction depends on the interaction range cutoff
:math:`l_{\mathrm{cut}}`, the Verlet list skin :math:`l_{\mathrm{skin}}`
and the box length :math:`l_{\mathrm{box}}`, and is determined automatically
by solving several equations. It can be useful to know how to estimate the
number of cells per direction, because it limits the number of MPI ranks
that can be allocated to an MPI-parallel simulation. As a rule of thumb,
for a cubic box geometry the number of cells per direction is often:

.. math::

    \left\lfloor \frac{l_{\mathrm{box}}}{l_{\mathrm{cut}} + l_{\mathrm{skin}}} \right\rfloor

For example, in a system with box length 12, LJ cutoff 2.5 and Verlet
skin 0.4, the number of cells cannot be more than 4 in each direction.
A runtime error will be triggered during integration when running a
simulation with such a system and allocating more than 64 MPI ranks
in total, or more than 4 MPI ranks per direction. In this situation,
consider increasing the box size or decreasing the interaction cutoff
or Verlet list skin.

.. _N-squared:

N-squared
^^^^^^^^^

Invoking :py:meth:`~espressomd.cell_system.CellSystem.set_n_square`
selects the very primitive N-squared cellsystem, which calculates
the interactions for all particle pairs. Therefore it loops over all
particles, giving an unfavorable computation time scaling of
:math:`N^2`. However, algorithms like MMM1D or the plain Coulomb
interaction in the cell model require the calculation of all pair
interactions. ::

    import espressomd
    system = espressomd.System(box_l=[1, 1, 1])
    system.cell_system.set_n_square()

In a multiple processor environment, the N-squared cellsystem uses a
simple particle balancing scheme to have a nearly equal number of
particles per CPU, :math:`n` nodes have :math:`m` particles, and
:math:`p-n` nodes have :math:`m+1` particles, such that
:math:`n \cdot m + (p - n) \cdot (m + 1) = N`, the total number of particles. Therefore the
computational load should be balanced fairly equal among the nodes, with
one exception: This code always uses one CPU for the interaction between
two different nodes. For an odd number of nodes, this is fine, because
the total number of interactions to calculate is a multiple of the
number of nodes, but for an even number of nodes, for each of the
:math:`p-1` communication rounds, one processor is idle.

E.g. for 2 processors, there are 3 interactions: 0-0, 1-1, 0-1.
Naturally, 0-0 and 1-1 are treated by processor 0 and 1, respectively.
But the 0-1 interaction is treated by node 1 alone, so the workload for
this node is twice as high. For 3 processors, the interactions are 0-0,
1-1, 2-2, 0-1, 1-2, 0-2. Of these interactions, node 0 treats 0-0 and
0-2, node 1 treats 1-1 and 0-1, and node 2 treats 2-2 and 1-2.

Therefore it is highly recommended that you use N-squared only with an
odd number of nodes, if with multiple processors at all.

.. _Hybrid decomposition:

Hybrid decomposition
^^^^^^^^^^^^^^^^^^^^

If for a simulation setup the interaction range is much smaller than the
system size, use of a :ref:`Regular decomposition` leads to efficient
scaling behavior (order :math:`N` instead of order :math:`N^2`).
Consider a system with many small particles, e.g. a polymer solution.
There, already the addition of one single large particle increases the maximum
interaction range and thus the minimum cell size of the decomposition.
Due to this larger cell size, throughout the simulation box a large number
of non-interacting pairs of small particles is visited during the short
range calculation. This can considerably increase the computational cost of
the simulation.

For such simulation setups, i.e. systems with a few large particles and much
more small particles, the hybrid decomposition can be used. This hybrid
decomposition is backed by two coupled particle decompositions which can
be used to efficiently deal with the differently sized particles.
Specifically that means putting the small particles into a
:ref:`Regular decomposition`. There, the minimum cell size is limited only
by the maximum interaction range of all particles within this decomposition.
The few large particles are put into a :ref:`N-squared` cellsystem. Particles
within this decomposition interact both, amongst each other and with all small
particles in the :ref:`Regular decomposition`. The hybrid decomposition can therefore
effectively recover the computational efficiency of the regular decomposition,
given that only a few large particles have been added.

Invoking :py:meth:`~espressomd.cell_system.CellSystem.set_hybrid_decomposition`
selects the hybrid decomposition. ::

    system = espressomd.System(box_l=[10, 10, 10])
    system.cell_system.set_hybrid_decomposition(n_square_types={1, 3}, cutoff_regular=1.2)

Here, ``n_square_types`` is a python set containing the types of particles to
put into the :ref:`N-squared` cellsystem, i.e. the particle types of the
large particles. Particles with other types will by default be put into the
:ref:`Regular decomposition`. Note that for now it is also necessary to manually set
the maximum cutoff to consider for interactions within the
:ref:`Regular decomposition`, i.e. the maximum interaction range among all
small particle types. Set this via the ``cutoff_regular`` parameter.

.. note::

  The hybrid particle decomposition has been added to |es| only recently and
  for now should be considered an experimental feature. If you notice some unexpected
  behavior please let us know via github or the mailing list.

