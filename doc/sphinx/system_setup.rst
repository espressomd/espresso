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

    * :py:attr:`~espressomd.cellsystem.CellSystem.node_grid`

    (int[3]) 3D node grid for real space domain decomposition (optional, if
    unset an optimal set is chosen automatically). The domain decomposition
    can be visualized with :file:`samples/visualization_cellsystem.py`.

    * :py:attr:`~espressomd.cellsystem.CellSystem.skin`

    (float) Skin for the Verlet list. This value has to be set, otherwise the simulation will not start.

Details about the cell system can be obtained by :meth:`espressomd.System().cell_system.get_state() <espressomd.cellsystem.CellSystem.get_state>`:

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

Invoking :py:meth:`~espressomd.cellsystem.CellSystem.set_n_square`
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

Therefore it is highly recommended that you use nsquared only with an
odd number of nodes, if with multiple processors at all.

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
:py:meth:`~espressomd.thermostat.Thermostat.set_langevin` of the thermostat
class :class:`espressomd.thermostat.Thermostat` has to be invoked.
Best explained in an example::

    import espressomd
    system = espressomd.System()
    therm = system.Thermostat()

    therm.set_langevin(kT=1.0, gamma=1.0, seed=41)

As explained before the temperature is set as thermal energy :math:`k_\mathrm{B} T`.

The Langevin thermostat is based on an extension of Newton's equation of motion to

.. math::  m_i \dot{v}_i(t) = f_i(\{x_j\},v_i,t) - \gamma v_i(t) + \sqrt{2\gamma k_B T} \eta_i(t).

Here, :math:`f_i` are all deterministic forces from interactions,
:math:`\gamma` the bare friction coefficient and :math:`\eta` a random, "thermal" force.
The friction term accounts for dissipation in a surrounding fluid whereas
the random force  mimics collisions of the particle with solvent molecules
at temperature :math:`T` and satisfies

.. math:: <\eta(t)> = 0 , <\eta^\alpha_i(t)\eta^\beta_j(t')> = \delta_{\alpha\beta} \delta_{ij}\delta(t-t')

(:math:`<\cdot>` denotes the ensemble average and :math:`\alpha,\beta` are spatial coordinates).

In the |es| implementation of the Langevin thermostat,
the additional terms only enter in the force calculation.
This reduces the accuracy of the Velocity Verlet integrator
by one order in :math:`dt` because forces are now velocity-dependent.

The random process :math:`\eta(t)` is discretized by drawing an uncorrelated random number
:math:`\overline{\eta}` for each component of all the particle forces.
The distribution of :math:`\overline{\eta}` is uniform and satisfies

.. math:: <\overline{\eta}> = 0 , <\overline{\eta}\overline{\eta}> = 1/dt

The keyword ``seed`` controls the state of the random number generator (Philox
Counter-based RNG) and is required on first activation of the thermostat. It
can be omitted in subsequent calls of ``set_langevin()``. It is the user's
responsibility to decide whether the thermostat should be deterministic (by
using a fixed seed) or not (by using a randomized seed).

If the feature ``ROTATION`` is compiled in, the rotational degrees of freedom are
also coupled to the thermostat. If only the first two arguments are
specified then the friction coefficient for the rotation is set to the
same value as that for the translation.
A separate rotational friction coefficient can be set by inputting
``gamma_rotate``. The two options allow one to switch the translational and rotational
thermalization on or off separately, maintaining the frictional behavior. This
can be useful, for instance, in high Péclet number active matter systems, where
one only wants to thermalize only the rotational degrees of freedom and
translational motion is effected by the self-propulsion.

The keywords ``gamma`` and ``gamma_rotate`` can be specified as a scalar,
or, with feature ``PARTICLE_ANISOTROPY`` compiled in, as the three eigenvalues
of the respective friction coefficient tensor. This is enables the simulation of
the anisotropic diffusion of anisotropic colloids (rods, etc.).

Using the Langevin thermostat, it is possible to set a temperature and a
friction coefficient for every particle individually via the feature
``LANGEVIN_PER_PARTICLE``.  Consult the reference of the ``part`` command
(chapter :ref:`Setting up particles`) for information on how to achieve this.

.. _LB thermostat:

Lattice-Boltzmann thermostat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :ref:`Lattice-Boltzmann` thermostat acts similar to the :ref:`Langevin thermostat` in that the governing equation for particles is

.. math::  m_i \dot{v}_i(t) = f_i(\{x_j\},v_i,t) - \gamma (v_i(t)-u(x_i(t),t)) + \sqrt{2\gamma k_B T} \eta_i(t).

where :math:`u(x,t)` is the fluid velocity at position :math:`x` and time :math:`t`.
To preserve momentum, an equal and opposite friction force and random force act on the fluid.

Numerically the fluid velocity is determined from the lattice-Boltzmann node velocities
by interpolating as described in :ref:`Interpolating velocities`.
The backcoupling of friction forces and noise to the fluid is also done by distributing those forces amongst the nearest LB nodes.
Details for both the interpolation and the force distribution can be found in :cite:`ahlrichs99` and :cite:`duenweg08a`.

The LB fluid can be used to thermalize particles, while also including their hydrodynamic interactions.
The LB thermostat expects an instance of either :class:`espressomd.lb.LBFluid` or :class:`espressomd.lb.LBFluidGPU`.
Temperature is set via the ``kT`` argument of the LB fluid.

Furthermore a ``seed`` has to be given for the
thermalization of the particle coupling. The magnitude of the frictional coupling can be adjusted by
the parameter ``gamma``.
To enable the LB thermostat, use::

    system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)


No other thermostatting mechanism is necessary
then. Please switch off any other thermostat before starting the LB
thermostatting mechanism.

The LBM implementation provides a fully thermalized LB fluid, all
nonconserved modes, including the pressure tensor, fluctuate correctly
according to the given temperature and the relaxation parameters. All
fluctuations can be switched off by setting the temperature to 0.

.. note:: Coupling between LB and MD only happens if the LB thermostat is set with a :math:`\gamma \ge 0.0`.


.. _Dissipative Particle Dynamics (DPD):

Dissipative Particle Dynamics (DPD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DPD thermostat adds friction and noise to the particle
dynamics like the :ref:`Langevin thermostat`, but these
are not applied to every particle individually but instead
encoded in a dissipative interaction between particles :cite:`soddeman03a`.

To realize a complete DPD fluid model in |es|, three parts are needed:
The DPD thermostat, which controls the temperate, a dissipative
interaction between the particles that make up the fluid,
see :ref:`DPD interaction`, and a repulsive conservative force.

The temperature is set via
:py:meth:`espressomd.thermostat.Thermostat.set_dpd`
which takes ``kT`` as the only argument.

The friction coefficients and cutoff are controlled via the
:ref:`DPD interaction` on a per type-pair basis. For details
see there.

The friction (dissipative) and noise (random) term are coupled via the
fluctuation-dissipation theorem. The friction term is a function of the
relative velocity of particle pairs. The DPD thermostat is better for
dynamics than the Langevin thermostat, since it mimics hydrodynamics in
the system.

As a conservative force any interaction potential can be used,
see :ref:`Isotropic non-bonded interactions`. A common choice is
a force ramp which is implemented as :ref:`Hat interaction`.

A complete example of setting up a DPD fluid and running it
to sample the equation of state can be found in samples/dpd.py.

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
Activate the NPT thermostat with the command :py:meth:`~espressomd.thermostat.Thermostat.set_npt`
and setup the integrator for the NPT ensemble with :py:meth:`~espressomd.integrate.IntegratorHandle.set_isotropic_npt`.

For example::

    import espressomd

    system = espressomd.System()
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0, seed=41)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

For an explanation of the algorithm involved, see :ref:`Isotropic NPT integrator`.

Be aware that this feature is neither properly examined for all systems
nor is it maintained regularly. If you use it and notice strange
behavior, please contribute to solving the problem.

.. _Brownian thermostat:

Brownian thermostat
~~~~~~~~~~~~~~~~~~~

Brownian thermostat is a formal name of a thermostat enabling the
Brownian Dynamics feature (see :cite:`schlick2010`) which implies
a propagation scheme involving systematic and thermal parts of the
classical Ermak-McCammom's (see :cite:`ermak78a`)
Brownian Dynamics. Currently it is implemented without
hydrodynamic interactions, i.e.
with a diagonal diffusion tensor.
The hydrodynamic interactions feature will be available later
as a part of the present Brownian Dynamics or
implemented separately within the Stokesian Dynamics.

In order to activate the Brownian thermostat, the member function
:py:attr:`~espressomd.thermostat.Thermostat.set_brownian` of the thermostat
class :class:`espressomd.thermostat.Thermostat` has to be invoked.
The system integrator should be also changed.
Best explained in an example::

    import espressomd
    system = espressomd.System()
    system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=41)
    system.integrator.set_brownian_dynamics()

where ``gamma`` (hereinafter :math:`\gamma`) is a viscous friction coefficient.
In terms of the Python interface and setup, the Brownian thermostat is very
similar to the :ref:`Langevin thermostat`. The feature
``BROWNIAN_PER_PARTICLE`` is used to control the per-particle
temperature and the friction coefficient setup. The major differences are
its internal integrator implementation and other temporal constraints.
The integrator is still a symplectic Velocity Verlet-like one.
It is implemented via a viscous drag part and a random walk of both the position and
velocity. Due to a nature of the Brownian Dynamics method, its time step :math:`\Delta t`
should be large enough compared to the relaxation time
:math:`m/\gamma` where :math:`m` is the particle mass.
This requirement is just a conceptual one
without specific implementation technical restrictions.
Note that with all similarities of
Langevin and Brownian Dynamics, the Langevin thermostat temporal constraint
is opposite. A velocity is restarting from zero at every step.
Formally, the previous step velocity at the beginning of the the :math:`\Delta t` interval
is dissipated further
and does not contribute to the end one as well as to the positional random walk.
Another temporal constraint
which is valid for both Langevin and Brownian Dynamics: conservative forces
should not change significantly over the :math:`\Delta t` interval.

The viscous terminal velocity :math:`\Delta v` and corresponding positional
step :math:`\Delta r` are fully driven by conservative forces :math:`F`:

.. math:: \Delta r = \frac{F \cdot \Delta t}{\gamma}

.. math:: \Delta v = \frac{F}{\gamma}

A positional random walk variance of each coordinate :math:`\sigma_p^2`
corresponds to a diffusion within the Wiener process:

.. math:: \sigma_p^2 = 2 \frac{kT}{\gamma} \cdot \Delta t

Each velocity component random walk variance :math:`\sigma_v^2` is defined by the heat
component:

.. math:: \sigma_v^2 = \frac{kT}{m}

Note that the velocity random walk is propagated from zero at each step.

A rotational motion is implemented similarly.
Note: the rotational Brownian dynamics implementation is compatible with particles which have
the isotropic moment of inertia tensor only. Otherwise, the viscous terminal angular velocity
is not defined, i.e. it has no constant direction over the time.

The keyword ``seed`` controls the state of the random number generator (Philox
Counter-based RNG) and is required on first activation of the thermostat. It
can be omitted in subsequent calls of ``set_brownian()``. It is the user's
responsibility to decide whether the thermostat should be deterministic (by
using a fixed seed) or not (by using a randomized seed).

.. _CUDA:

CUDA
----

:py:meth:`~espressomd.cuda_init.CudaInitHandle()` command can be used to choose the GPU for all subsequent
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
the device id as key and the device name as its value.

.. _Selection of CUDA device:

Selection of CUDA device
~~~~~~~~~~~~~~~~~~~~~~~~

When you start ``pypresso`` your first GPU should be selected.
If you wanted to use the second GPU, this can be done
by setting :attr:`espressomd.cuda_init.CudaInitHandle.device` as follows::

    system = espressomd.System()

    system.cuda_init_handle.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.
