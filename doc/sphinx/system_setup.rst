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
    :py:func:`~espressomd.system.System.change_volume_and_rescale_particles`.

* :py:attr:`~espressomd.system.System.periodicity`

    (int[3]) Specifies periodicity for the three directions. If the feature
    ``PARTIAL_PERIODIC`` is set, |es| can be instructed to treat some
    dimensions as non-periodic. Per default espresso assumes periodicity in
    all directions which equals setting this variable to [1,1,1]. A
    dimension is specified as non-periodic via setting the periodicity
    variable for this dimension to 0. E.g. Periodicity only in z-direction
    is obtained by [0,0,1]. Caveat: Be aware of the fact that making a
    dimension non-periodic does not hinder particles from leaving the box in
    this direction. In this case for keeping particles in the simulation box
    a constraint has to be set.

* :py:attr:`~espressomd.system.System.time_step`

    (float) Time step for MD integration.

* :py:attr:`~espressomd.system.System.time`

    (float) The simulation time.

* :py:attr:`~espressomd.system.System.min_global_cut`

    (float) Minimal total cutoff for real space. Effectively, this plus the
    :py:attr:`~espressomd.cellsystem.CellSystem.skin` is the minimally possible cell size. Espresso typically determines
    this value automatically, but some algorithms, virtual sites, require
    you to specify it manually.

* :py:attr:`~espressomd.system.System.max_cut_bonded`

    *read-only* Maximal cutoff of bonded real space interactions.

* :py:attr:`~espressomd.system.System.max_cut_nonbonded`

    *read-only* Maximal cutoff of bonded real space interactions.

.. _Accessing module states:

Accessing module states
~~~~~~~~~~~~~~~~~~~~~~~

Some variables like or are no longer directly available as attributes.
In these cases they can be easily derived from the corresponding Python
objects like

``n_part = len(espressomd.System().part[:].pos)``

or by calling the corresponding ``get_state()`` methods like::

    temperature = espressomd.System().thermostat.get_state()[0]['kT']

    gamma = espressomd.System().thermostat.get_state()[0]['gamma']

    gamma_rot = espressomd.System().thermostat.get_state()[0]['gamma_rotation']

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

    * :py:attr:`~espressomd.cellsystem.CellSystem.max_num_cells`

    (int) Maximal number of cells for the link cell algorithm. Reasonable
    values are between 125 and 1000, or for some problems :math:`n_part / nnodes`.

    * :py:attr:`~espressomd.cellsystem.CellSystem.min_num_cells`

    (int) Minimal number of cells for the link cell algorithm. Reasonable
    values range in :math:`10^{-6} N^2` to :math:`10^{-7} N^2`. In general
    just make sure that the Verlet lists are not incredibly large. By default the
    minimum is 0, but for the automatic P3M tuning it may be wise to set larger
    values for high particle numbers.

    * :py:attr:`~espressomd.cellsystem.CellSystem.node_grid`

    (int[3]) 3D node grid for real space domain decomposition (optional, if
    unset an optimal set is chosen automatically).

    * :py:attr:`~espressomd.cellsystem.CellSystem.skin`

    (float) Skin for the Verlet list. This value has to be set, otherwise the simulation will not start.

Details about the cell system can be obtained by :meth:`espressomd.System().cell_system.get_state() <espressomd.cellsystem.CellSystem.get_state>`:

    * ``cell_grid``       Dimension of the inner cell grid.
    * ``cell_size``       Box-length of a cell.
    * ``local_box_l``     Local simulation box length of the nodes.
    * ``max_cut``         Maximal cutoff of real space interactions.
    * ``n_layers``        Number of layers in cell structure LAYERED
    * ``n_nodes``         Number of nodes.
    * ``type``            The current type of the cell system.
    * ``verlet_reuse``    Average number of integration steps the Verlet list is re-used.

.. _Domain decomposition:

Domain decomposition
~~~~~~~~~~~~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_domain_decomposition`
selects the domain decomposition cell scheme, using Verlet lists
for the calculation of the interactions. If you specify ``use_verlet_lists=False``, only the
domain decomposition is used, but not the Verlet lists. ::

    system = espressomd.System()

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
therefore of the order N instead of order :math:`N^2` if one has to
calculate all pair interactions.

.. _N-squared:

N-squared
~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_n_square`
selects the very primitive nsquared cellsystem, which calculates
the interactions for all particle pairs. Therefore it loops over all
particles, giving an unfavorable computation time scaling of
:math:`N^2`. However, algorithms like MMM1D or the plain Coulomb
interaction in the cell model require the calculation of all pair
interactions. ::

    system = espressomd.System()

    system.cell_system.set_n_square()

In a multiple processor environment, the nsquared cellsystem uses a
simple particle balancing scheme to have a nearly equal number of
particles per CPU, :math:`n` nodes have :math:`m` particles, and
:math:`p-n` nodes have :math:`m+1` particles, such that
:math:`n*m+(p-n)*(m+1)=N`, the total number of particles. Therefore the
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

Therefore it is highly recommended that you use nsquared only with an
odd number of nodes, if with multiple processors at all.

.. _Layered cell system:

Layered cell system
~~~~~~~~~~~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_layered`
selects the layered cell system, which is specifically designed for
the needs of the MMM2D algorithm. Basically it consists of a nsquared
algorithm in x and y, but a domain decomposition along z, i.e. the
system is cut into equally sized layers along the z axis. The current
implementation allows for the CPUs to align only along the z axis,
therefore the processor grid has to have the form 1x1xN. However, each
processor may be responsible for several layers, which is determined by
``n_layers``, i.e. the system is split into N\* layers along the z axis. Since in x
and y direction there are no processor boundaries, the implementation is
basically just a stripped down version of the domain decomposition
cellsystem. ::

    system = espressomd.System()

    system.cell_system.set_layered(n_layers=4)

.. _Thermostats:

Thermostats
-----------

The thermostat can be controlled by the class :class:`espressomd.thermostat.Thermostat`.

The different available thermostats will be described in the following
subsections. Note that for a simulation of the NPT ensemble, you need to
use a standard thermostat for the particle velocities (Langevin or DPD),
and a thermostat for the box geometry (the isotropic NPT thermostat).

You may combine different thermostats at your own risk by turning them
on one by one. Note that there is only one temperature for all
thermostats, although for some thermostats like the Langevin thermostat,
particles can be assigned individual temperatures.

Since |es| does not enforce a particular unit system, it cannot know about
the current value of the Boltzmann constant. Therefore, when specifying
the temperature of a thermostat, you actually do not define the
temperature, but the value of the thermal energy :math:`k_B T` in the
current unit system (see the discussion on units, Section :ref:`On units`).

.. _Langevin thermostat:

Langevin thermostat
~~~~~~~~~~~~~~~~~~~

In order to activate the Langevin thermostat the member function
:py:attr:`~espressomd.thermostat.Thermostat.set_langevin` of the thermostat
class :class:`espressomd.thermostat.Thermostat` has to be invoked.
Best explained in an example::

    import espressomd
    system = espressomd.System()
    therm = system.Thermostat()

    therm.set_langevin(kT=1.0, gamma=1.0, seed=41)

As explained before the temperature is set as thermal energy :math:`k_\mathrm{B} T`.
The Langevin thermostat consists of a friction and noise term coupled
via the fluctuation-dissipation theorem. The friction term is a function
of the particle velocities. By specifying the diffusion coefficient for
the particle becomes

.. math:: D = \frac{\text{temperature}}{\text{gamma}}.

The relaxation time is given by :math:`\text{gamma}/\text{MASS}`, with
``MASS`` the particle's mass.  For a more detailed explanation, refer to
:cite:`grest86a`.  An anisotropic diffusion coefficient tensor is available to
simulate anisotropic colloids (rods, etc.) properly. It can be enabled by the
feature ``PARTICLE_ANISOTROPY``.

The keyword ``seed`` controls the state of the random number generator (Philox
Counter-based RNG) and is required on first activation of the thermostat. It
can be omitted in subsequent calls of ``set_langevin()``. It is the user's
responsibility to decide whether the thermostat should be deterministic (by
using a fixed seed) or not (by using a randomized seed).

If the feature ``ROTATION`` is compiled in, the rotational degrees of freedom are
also coupled to the thermostat. If only the first two arguments are
specified then the diffusion coefficient for the rotation is set to the
same value as that for the translation.

A separate rotational diffusion coefficient can be set by inputting
``gamma_rotate``.  This also allows one to properly match the translational and
rotational diffusion coefficients of a sphere. Feature ``ROTATIONAL_INERTIA``
enables an anisotropic rotational diffusion coefficient tensor through
corresponding friction coefficients.

Finally, the two options allow one to switch the translational and rotational
thermalization on or off separately, maintaining the frictional behavior. This
can be useful, for instance, in high Péclet number active matter systems, where
one only wants to thermalize the rotational degrees of freedom and
translational motion is effected by the self-propulsion.

Using the Langevin thermostat, it is possible to set a temperature and a
friction coefficient for every particle individually via the feature
``LANGEVIN_PER_PARTICLE``.  Consult the reference of the ``part`` command
(chapter :ref:`Setting up particles`) for information on how to achieve this.

.. _Dissipative Particle Dynamics (DPD):

Dissipative Particle Dynamics (DPD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To realize a complete DPD fluid model, three parts are needed:
The DPD thermostat, which controls the temperate, a dissipative
interaction between the particles that make up the fluid,
see :ref:`DPD interaction`, and a repulsive conservative force.

The DPD thermostat can be invoked by the function:
:py:attr:`espressomd.thermostat.Thermostat.set_dpd`
which takes :math:`k_\mathrm{B} T` as the only argument.

The friction coefficients and cutoff are controlled via the
:ref:`DPD interaction` on a per type-pair basis. For details see
there.

As a conservative force any interaction potential can be used,
see :ref:`Isotropic non-bonded interactions`. A common choice is
a force ramp which is implemented as :ref:`Hat interaction`.

A complete example of setting up a DPD fluid and running it
to sample the equation of state can be found in samples/dpd.py.

DPD adds a velocity dependent dissipative force and a random force
to the conservative pair forces.

The friction (dissipative) and noise (random) term are coupled via the
fluctuation- dissipation theorem. The friction term is a function of the
relative velocity of particle pairs. The DPD thermostat is better for
dynamics than the Langevin thermostat, since it mimics hydrodynamics in
the system.

When using a Lennard-Jones interaction, :math:`{r_\mathrm{cut}} =
2^{\frac{1}{6}} \sigma` is a good value to choose, so that the
thermostat acts on the relative velocities between nearest neighbor
particles. Larger cutoffs including next nearest neighbors or even more
are unphysical.

Boundary conditions for DPD can be introduced by adding the boundary
as a particle constraint, and setting a velocity and a type on it, see
:class:`espressomd.constraints.Constraint`. Then a
:ref:`DPD interaction` with the type can be defined, which acts as a
boundary condition.

.. _Isotropic NPT thermostat:

Isotropic NPT thermostat
~~~~~~~~~~~~~~~~~~~~~~~~

This feature allows to simulate an (on average) homogeneous and isotropic system in the NPT ensemble.
In order to use this feature, ``NPT`` has to be defined in the :file:`myconfig.hpp`.
Activate the NPT thermostat with the command :py:func:`~espressomd.thermostat.Thermostat.set_npt`
and set the following parameters:

    * ``kT``:     (float) Thermal energy of the heat bath
    * ``gamma0``: (float) Friction coefficient of the bath
    * ``gammav``: (float) Artificial friction coefficient for the volume fluctuations.

Also, setup the integrator for the NPT ensemble with :py:func:`~espressomd.system.integrator.set_isotropic_npt`
and the parameters:

    * ``ext_pressure``:  (float) The external pressure as float variable.
    * ``piston``:        (float) The mass of the applied piston as float variable.

This thermostat is based on the Andersen thermostat (see
:cite:`andersen80a,mann05d`) and will thermalize the box
geometry. It will only do isotropic changes of the box.
See this code snippet for the two commands::

    import espressomd

    system = espressomd.System()
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

Be aware that this feature is neither properly examined for all systems
nor is it maintained regularly. If you use it and notice strange
behavior, please contribute to solving the problem.

.. _CUDA:

CUDA
----

:py:attr:`~espressomd.cuda_init.CudaInitHandle()` command can be used to choose the GPU for all subsequent
GPU-computations. Note that due to driver limitations, the GPU cannot be
changed anymore after the first GPU-using command has been issued, for
example ``lbfluid``. If you do not choose the GPU manually before that,
CUDA internally chooses one, which is normally the most powerful GPU
available, but load-independent. ::

    system = espressomd.System()

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

If you want to list available CUDA devices
you should access :attr:`espressomd.cuda_init.CudaInitHandle.device_list`, e.g., ::

    system = espressomd.System()

    print(system.cuda_init_handle.device_list)

This attribute is read only and will return a dictionary containing
the device id as key and the device name as its' value.

.. _Selection of CUDA device:

Selection of CUDA device
~~~~~~~~~~~~~~~~~~~~~~~~

When you start ``pypresso`` your first GPU should
be selected.
If you wanted to use the second GPU, this can be done
by setting :attr:`espressomd.cuda_init.CudaInitHandle.device` as follows::

    system = espressomd.System()

    system.cuda_init_handle.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.
