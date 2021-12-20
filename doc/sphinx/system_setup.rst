.. _Setting up the system:

Setting up the system
=====================

.. _Setting global variables in Python:

Setting global variables in Python
----------------------------------

The global variables in Python are controlled via the
:class:`espressomd.system.System` class.
Global system variables can be read and set in Python simply by accessing the
attribute of the corresponding Python object. Those variables that are already
available in the Python interface are listed in the following. Note that for the
vectorial properties ``box_l`` and ``periodicity``, component-wise manipulation
like ``system.box_l[0] = 1`` or in-place operators like ``+=`` or ``*=`` are not
allowed and result in an error. This behavior is inherited, so the same applies
to ``a`` after ``a = system.box_l``. If you want to use a vectorial property
for further calculations, you should explicitly make a copy e.g. via
``a = numpy.copy(system.box_l)``.

* :py:attr:`~espressomd.system.System.box_l`

  (float[3]) Simulation box lengths of the cuboid box used by |es|.
  Note that if you change the box length during the simulation, the folded
  particle coordinates will remain the same, i.e., the particle stay in
  the same image box, but at the same relative position in their image
  box. If you want to scale the positions, use the command
  :py:meth:`~espressomd.system.System.change_volume_and_rescale_particles`.

* :py:attr:`~espressomd.system.System.periodicity`

  (int[3]) Specifies periodicity for the three directions. |es| can be instructed
  to treat some dimensions as non-periodic. By default |es| assumes periodicity in
  all directions which equals setting this variable to ``[True, True, True]``.
  A dimension is specified as non-periodic via setting the periodicity
  variable for this dimension to ``False``. E.g. Periodicity only in z-direction
  is obtained by ``[False, False, True]``. Caveat: Be aware of the fact that making a
  dimension non-periodic does not hinder particles from leaving the box in
  this direction. In this case for keeping particles in the simulation box
  a constraint has to be set.

* :py:attr:`~espressomd.system.System.time_step`

  (float) Time step for MD integration.

* :py:attr:`~espressomd.system.System.time`

  (float) The simulation time.

* :py:attr:`~espressomd.system.System.min_global_cut`

  (float) Minimal total cutoff for real space. Effectively, this plus the
  :py:attr:`~espressomd.cellsystem.CellSystem.skin` is the minimally possible
  cell size. |es| typically determines this value automatically, but some
  algorithms, virtual sites, require you to specify it manually.

* :py:attr:`~espressomd.system.System.max_cut_bonded`

  *read-only* Maximal cutoff of bonded real space interactions.

* :py:attr:`~espressomd.system.System.max_cut_nonbonded`

  *read-only* Maximal cutoff of bonded real space interactions.

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

.. _Cellsystems:

Cellsystems
-----------

This section deals with the flexible particle data organization of |es|. Due
to different needs of different algorithms, |es| is able to change the
organization of the particles in the computer memory, according to the
needs of the used algorithms. For details on the internal organization,
refer to section :ref:`Internal particle organization`.

.. _Global properties:

Global properties
~~~~~~~~~~~~~~~~~

The properties of the cell system can be accessed by
:class:`espressomd.system.System.cell_system`:

* :py:attr:`~espressomd.cellsystem.CellSystem.node_grid`

  (int[3]) 3D node grid for real space domain decomposition (optional, if
  unset an optimal set is chosen automatically). The domain decomposition
  can be visualized with :file:`samples/visualization_cellsystem.py`.

* :py:attr:`~espressomd.cellsystem.CellSystem.skin`

  (float) Skin for the Verlet list. This value has to be set, otherwise the simulation will not start.

Details about the cell system can be obtained by :meth:`espressomd.system.System.cell_system.get_state() <espressomd.cellsystem.CellSystem.get_state>`:

* ``cell_grid``       Dimension of the inner cell grid.
* ``cell_size``       Box-length of a cell.
* ``local_box_l``     Local simulation box length of the nodes.
* ``max_cut``         Maximal cutoff of real space interactions.
* ``n_nodes``         Number of nodes.
* ``type``            The current type of the cell system.
* ``verlet_reuse``    Average number of integration steps the Verlet list is re-used.

.. _Domain decomposition:

Domain decomposition
~~~~~~~~~~~~~~~~~~~~

Invoking :py:meth:`~espressomd.cellsystem.CellSystem.set_domain_decomposition`
selects the domain decomposition cell scheme, using Verlet lists
for the calculation of the interactions. If you specify ``use_verlet_lists=False``, only the
domain decomposition is used, but not the Verlet lists. ::

    system = espressomd.System(box_l=[1, 1, 1])

    system.cell_system.set_domain_decomposition(use_verlet_lists=True)

The domain decomposition cellsystem is the default system and suits most
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

.. _N-squared:

N-squared
~~~~~~~~~

Invoking :py:meth:`~espressomd.cellsystem.CellSystem.set_n_square`
selects the very primitive N-squared cellsystem, which calculates
the interactions for all particle pairs. Therefore it loops over all
particles, giving an unfavorable computation time scaling of
:math:`N^2`. However, algorithms like MMM1D or the plain Coulomb
interaction in the cell model require the calculation of all pair
interactions. ::

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


.. _CUDA:

CUDA
----

:py:meth:`~espressomd.cuda_init.CudaInitHandle()` command can be used to choose the GPU for all subsequent
GPU-computations. Note that due to driver limitations, the GPU cannot be
changed anymore after the first GPU-using command has been issued, for
example ``lbfluid``. If you do not choose the GPU manually before that,
CUDA internally chooses one, which is normally the most powerful GPU
available, but load-independent. ::

    system = espressomd.System(box_l=[1, 1, 1])
    dev = system.cuda_init_handle.device
    system.cuda_init_handle.device = dev

The first invocation in the sample above returns the id of the set graphics card, the second one sets the
device id.

.. _GPU Acceleration with CUDA:

GPU Acceleration with CUDA
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``CUDA`` required

|es| is capable of GPU acceleration to speed up simulations.
Not every simulation method is parallelizable or profits from
GPU acceleration. Refer to :ref:`Available simulation methods`
to check whether your desired method can be used on the GPU.
In order to use GPU acceleration you need a NVIDIA GPU
and it needs to have at least compute capability 2.0.

For more information please check :class:`espressomd.cuda_init.CudaInitHandle`.

.. _List available CUDA devices:

List available CUDA devices
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to list available CUDA devices, you should call
:meth:`espressomd.cuda_init.CudaInitHandle.list_devices`::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> print(system.cuda_init_handle.list_devices())
    {0: 'GeForce RTX 2080', 1: 'GeForce GT 730'}

This method returns a dictionary containing
the device id as key and the device name as its value.

To get more details on the CUDA devices for each MPI node, call
:meth:`espressomd.cuda_init.CudaInitHandle.list_devices_properties`::

    >>> import pprint
    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> pprint.pprint(system.cuda_init_handle.list_devices_properties())
    {'seraue': {0: {'name': 'GeForce RTX 2080',
                    'compute_capability': (7, 5),
                    'cores': 46,
                    'total_memory': 8370061312},
                1: {'name': 'GeForce GT 730',
                    'compute_capability': (3, 5),
                    'cores': 2,
                    'total_memory': 1014104064}}}

.. _Selection of CUDA device:

Selection of CUDA device
~~~~~~~~~~~~~~~~~~~~~~~~

When you start ``pypresso`` your first GPU should be selected.
If you wanted to use the second GPU, this can be done
by setting :attr:`espressomd.cuda_init.CudaInitHandle.device` as follows::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> system.cuda_init_handle.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.
