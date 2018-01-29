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
to ``a`` after ``a = system.box_l``. If you want to use an vectorial property
for further calculations, you should explicity make a copy e.g. via
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
    PARTIAL\_PERIODIC is set, |es| can be instructed to treat some
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


Accessing module states
~~~~~~~~~~~~~~~~~~~~~~~

Some variables like or are no longer directly available as attributes.
In these cases they can be easily derived from the corresponding Python
objects like

``n_part = len(espressomd.System().part[:].pos)``

or by calling the corresponding ``get_state`` methods like::

    temperature = espressomd.System().thermostat.get_state()[0][’kT’]
    
    gamma = espressomd.System().thermostat.get_state()[0][’gamma’]
    
    gamma_rot = espressomd.System().thermostat.get_state()[0][’gamma_rotation’]

.. _\`\`thermostat\`\`\: Setting up the thermostat:

``thermostat``: Setting up the thermostat
-----------------------------------------

The thermostat can be controlled by the class :class:`espressomd.thermostat.Thermostat`

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
current unit system (see the discussion on units, Section [sec:units]).

Note that there are three different types of noise which can be used in
|es|. The one used typically in simulations is flat noise with the correct
variance and it is the default used in |es|, though it can be explicitly
specified using the feature ``FLATNOISE``. You can also employ Gaussian noise which
is, in some sense, more realistic. Notably Gaussian noise (activated
using the feature ``GAUSSRANDOM``) does a superior job of reproducing higher order
moments of the Maxwell-Boltzmann distribution. For typical generic
coarse-grained polymers using FENE bonds the Gaussian noise tends to
break the FENE bonds. We thus offer a third type of noise, activate
using the feature ``GAUSSRANDOMCUT``, which produces Gaussian random numbers but takes
anything which is two standard deviations (:math:`2\sigma`) below or
above zero and set it to :math:`-2\sigma` or :math:`2\sigma`
respectively. In all three cases the distribution is made such that the
second moment of the distribution is the same and thus results in the
same temperature.

.. _Langevin thermostat:

Langevin thermostat
~~~~~~~~~~~~~~~~~~~

In order to activate the Langevin thermostat the member function
:py:attr:`~espressomd.thermostat.Thermostat.set_langevin` of the thermostat
class :class:`espressomd.thermostat.Thermostat` has to be invoked.
Best explained in an example:::
    
    import espressomd
    system = espressomd.System()
    therm  = system.Thermostat()

    therm.set_langevin(kT=1.0, gamma=1.0)

As explained before the temperature is set as thermal energy :math:`k_\mathrm{B} T`. 
The Langevin thermostat consists of a friction and noise term coupled
via the fluctuation-dissipation theorem. The friction term is a function
of the particle velocities. By specifying the diffusion coefficient for
the particle becomes

.. math:: D = \frac{\text{temperature}}{\text{gamma}}.

The relaxation time is given by :math:`\text{gamma}/\text{MASS}`, with
``MASS`` the particle’s mass.  For a more detailed explanation, refer to
:cite:`grest86a`.  An anisotropic diffusion coefficient tensor is available to
simulate anisotropic colloids (rods, etc.) properly. It can be enabled by the
feature ``PARTICLE_ANISOTROPY``.

If the feature ``ROTATION`` is compiled in, the rotational degrees of freedom are
also coupled to the thermostat. If only the first two arguments are
specified then the diffusion coefficient for the rotation is set to the
same value as that for the translation.

A separate rotational diffusion coefficient can be set by inputting
``gamma_rotate``.  This also allows one to properly match the translational and
rotational diffusion coefficients of a sphere. ``ROTATIONAL_INERTIA`` Feature
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

GHMC thermostat
~~~~~~~~~~~~~~~

.. todo::
    this is not yet implemented in the python interface.


Implements Generalized Hybrid Monte Carlo (GHMC) as a thermostat. GHMC
is a simulation method for sampling the canonical ensemble
:cite:`mehlig92`. The method consists of MC cycles that
combine a few constant energy MD steps, specified by , followed by a
Metropolis criterion for their acceptance. Prior to integration, the
particles momenta are mixed with momenta sampled from the appropriate
Boltzmann distribution.

Given the particles momenta :math:`\mathbf{p}^j` from the last
:math:`j^{th}` GHMC cycle the new momenta are generated by:
:math:`\mathbf{p}^{j+1}=\cos(\phi)\mathbf{p}^j+\sin(\phi)\pmb{\xi}`,
where :math:`\pmb{\xi}` is a noise vector of random Gaussian variables
with zero mean and variance :math:`1/\mathrm{temperature}` (see
:cite:`horowitz91` for more details). The momenta mixing
parameter :math:`\cos(\phi)` corresponds to in the implementation.

In case the MD step is rejected, the particles momenta may be flipped.
This is specified by setting the / option, for the option half of the
rejected MD steps randomly result in momenta flip. The default for
momenta flip is . The :math:`\pmb{\xi}` noise vector’s variance van be
tuned to exactly :math:`1/\mathrm{temperature}` by specifying the option.
The default for temperature scaling is .

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

When using a Lennard-Jones interaction, :math:`{r\_cut} =
2^{\frac{1}{6}} \sigma` is a good value to choose, so that the
thermostat acts on the relative velocities between nearest neighbor
particles. Larger cutoffs including next nearest neighbors or even more
are unphysical.

Isotropic NPT thermostat
~~~~~~~~~~~~~~~~~~~~~~~~

In order to use this feature, ``NPT`` has to be defined in the ``myconfig.hpp``.
Activate the NPT thermostat with the command :py:func:`~espressomd.thermostat.Thermostat.set_npt`
and set the following parameters:

    * kT:     (float) Thermal energy of the heat bath
    * gamma0: (float) Friction coefficient of the bath
    * gammav: (float) Artificial friction coefficient for the volume fluctuations.

Also, setup the integrator for the NPT ensemble with :py:func:`~espressomd.system.integrator.set_isotropic_npt` 
and the parameters:

    * ext_pressure:  (float) The external pressure as float variable.
    * piston:        (float) The mass of the applied piston as float variable.

This thermostat is based on the Anderson thermostat (see
:cite:`andersen80a,mann05d`) and will thermalize the box
geometry. It will only do isotropic changes of the box. 
See this code snippet for the two commands::

    import espressomd

    system=espressomd.System()
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

Be aware that this feature is neither properly examined for all systems
nor is it maintained regularly. If you use it and notice strange
behaviour, please contribute to solving the problem.

.. _\`\`nemd\`\`\: Setting up non-equilibirum MD:

``nemd``: Setting up non-equilibrium MD
---------------------------------------

.. todo::
    This is not implemented for the python interface yet

nemd exchange nemd shearrate nemd off nemd nemd profile nemd viscosity

Use NEMD (Non Equilibrium Molecular Dynamics) to simulate a system under
shear with help of an unphysical momentum change in two slabs in the
system.

Variants and will initialise NEMD. Two distinct methods exist. Both
methods divide the simulation box into slabs that lie parallel to the
x-y-plane and apply a shear in x direction. The shear is applied in the
top and the middle slabs. Note, that the methods should be used with a
DPD thermostat or in an NVE ensemble. Furthermore, you should not use
other special features like or inside the top and middle slabs. For
further reference on how NEMD is implemented into see
:cite:`soddeman01a`.

Variant chooses the momentum exchange method. In this method, in each
step the largest positive x-components of the velocity in the middle
slab are selected and exchanged with the largest negative x-components
of the velocity in the top slab.

Variant chooses the shear-rate method. In this method, the targeted
x-component of the mean velocity in the top and middle slabs are given
by

.. math:: {target\_velocity} = \pm {shearrate}\,\frac{L_z}{4}

where :math:`L_z` is the simulation box size in z-direction. During the
integration, the x-component of the mean velocities of the top and
middle slabs are measured. Then, the difference between the mean
x-velocities and the target x-velocities are added to the x-component of
the velocities of the particles in the respective slabs.

Variant will turn off NEMD, variant will print usage information of the
parameters of NEMD. Variant will return the velocity profile of the
system in x-direction (mean velocity per slab).

Variant will return the viscosity of the system, that is computed via

.. math:: \eta = \frac{F}{\dot{\gamma} L_x L_y}

where :math:`F` is the mean force (momentum transfer per unit time)
acting on the slab, :math:`L_x L_y` is the area of the slab and
:math:`\dot{\gamma}` is the shearrate.

NEMD as implemented generates a Poiseuille flow, with shear flow rate
varying over a finite wavelength determined by the box. For a planar
Couette flow (constant shear, infinite wavelength), consider using
Lees-Edwards boundary conditions (see ) to drive the shear.

.. _Lees-Edwards boundary conditions:

Lees-Edwards boundary condition: Setting up a shear flow
-------------------------------------------------------------

To use the Lees-Edwards boundary conditions, the feature ``LEES_EDWARDS`` is required.

Lees-Edwards boundary conditions can be used to introduce a shear flow to the MD simulation. An introduction can be found in :cite:`lees72`. Compared to NEMD simulations they have two big advantages: First, the bulk behavior of the system remains unchanged. Second, the image boxes are moved, whereas the flow within the primary simulation box has to develop on its own. Hence, this allows two additional phenomena: Shear banding can occur as well as non-linear shear profiles can be observed. This makes Lees-Edwards boundary conditions suitable for comparison with rheological experiments. 

Lees-Edwards boundary conditions impose a shear flow of speed :math:`\dot\gamma` by moving the periodic image boxes along the x-direction according to:

.. math:: v_{\text{x, unfolded}} = v_{\text{x, folded}} + \dot\gamma \cdot y_{\text{imagecount}}

:math:`v_{\text{x, unfolded}}` refers to the velocity of a particle outside the main simulation box, :math:`y_{\text{imagecount}}` is the amount of periodic boundaries crossed in the  :math:`y`-direction. 

The absolute offset of the periodic images can be set via

* :py:attr:`~espressomd.System().lees_edwards_offset`

The following example introduces the usage::
    
    import espressomd
    system = espressomd.System()
    absolute_offset = 0.2
    system.lees_edwards_offset = absolute_offset

Lees-Edwards boundary conditions can be used to obtain the shear modulus :math:`G = \frac{\tau}{\gamma}` or the shear viscosity :math:`\eta = \frac{\tau}{\dot\gamma}` outside the linear regime, where Green-Kubo relations are not valid anymore. For this purpose a lees_edwards_offset is set followed by one integration step for multiple times. Strain, strain rate and the shear stress need to be recorded for the calculation. Alternatively a sinusoidal lees_edwards_offset series can be used to carry out oscillatory experiments to calculate viscoelastic moduli (:math:`G', G''`). Furthermore a lees_edwards_offset can be set followed by many integration steps obtain the relaxation behaviour of a system. 

When applying a constant shear rate :math:`\dot\gamma` the velocity of the particles changes from :math:`-\frac{\dot\gamma}{2}` at the bottom of the box to :math:`\frac{\dot\gamma}{2}` at the top of the box. 

Physical meaningful values for systems where hydrodynamics play a major role, can only be obtained by including hydrodynamic interactions. Lees-Edwards boundary conditions are implemented in the :ref:`Lattice-Boltzmann` algorithms. For this algorithm the feature ``LB_GPU`` is required. Please refer to chapter :ref:`Lattice-Boltzmann` for more information. 

Lees-Edwards boundary conditions work with the DPD thermostat. In order to correctly observe transport properties, symmetry-breaking or entropy production in relation to shear flow is probably better to use the DPD thermostat (:ref:`Dissipative Particle Dynamics (DPD)`) once the initial heat-up has been carried out. The DPD thermostat removes kinetic energy from the system based on a frictional term defined relative to a local reference frame of a given particle-pair, without enforcing any specific flow pattern apriori. At high rates of dissipation, this can however lead to an artefactual shear-banding type effect at the periodic boundaries, such that the bulk fluid is nearly stationary. y. This effect is removed using the modification proposed to the DPD thermostat by Chatterjee :cite:`chatterjee2007` to allow treatment of systems with high dissipation rates, which is applied automatically if ``LEES_EDWARDS`` is compiled in. Chatterjee’s modification is just to skip calculation of DPD forces (both dissipative and random) for particle pairs which cross a boundary in y.

The command::

  print(system.lees_edwards_offset)

returns the current value of the offset. If ``LEES_EDWARDS`` is compiled in, then coordinates are folded into the primary simulation box as the integration progresses, to prevent a numerical overflow.

.. _cellsystem:

``cellsystem``: Setting up the cell system
------------------------------------------

This section deals with the flexible particle data organization of |es|. Due
to different needs of different algorithms, |es| is able to change the
organization of the particles in the computer memory, according to the
needs of the used algorithms. For details on the internal organization,
refer to section :ref:`internal_particle_org`.

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

Details about the cell system can be obtained by ``espressomd.System().cell_system.get_state()``:

    * `cell_grid`       Dimension of the inner cell grid.
    * `cell_size`       Box-length of a cell.
    * `local_box_l`     Local simulation box length of the nodes.
    * `max_cut`         Maximal cutoff of real space interactions.
    * `n_layers`        Number of layers in cell structure LAYERED
    * `n_nodes`         Number of nodes.
    * `type`            The current type of the cell system.
    * `verlet_reuse`    Average number of integration steps the verlet list is re-used.


Domain decomposition
~~~~~~~~~~~~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_domain_decomposition` 
selects the domain decomposition cell scheme, using Verlet lists
for the calculation of the interactions. If you specify ``use_verlet_lists=False``, only the
domain decomposition is used, but not the Verlet lists.::

    system=espressomd.System()

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

N-squared
~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_n_square`
selects the very primitive nsquared cellsystem, which calculates
the interactions for all particle pairs. Therefore it loops over all
particles, giving an unfavorable computation time scaling of
:math:`N^2`. However, algorithms like MMM1D or the plain Coulomb
interaction in the cell model require the calculation of all pair
interactions.::

    system=espressomd.System()

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

Layered cell system
~~~~~~~~~~~~~~~~~~~

Invoking :py:attr:`~espressomd.cellsystem.CellSystem.set_layered`
selects the layered cell system, which is specifically designed for
the needs of the MMM2D algorithm. Basically it consists of a nsquared
algorithm in x and y, but a domain decomposition along z, i. e. the
system is cut into equally sized layers along the z axis. The current
implementation allows for the CPUs to align only along the z axis,
therefore the processor grid has to have the form 1x1xN. However, each
processor may be responsible for several layers, which is determined by
``n_layers``, i. e. the system is split into N\* layers along the z axis. Since in x
and y direction there are no processor boundaries, the implementation is
basically just a stripped down version of the domain decomposition
cellsystem.::

    system=espressomd.System()

    system.cell_system.set_layered(n_layers=4)

CUDA
----

:py:attr:`~espressomd.cuda_init.CudaInitHandle()` command can be used to choose the GPU for all subsequent
GPU-computations. Note that due to driver limitations, the GPU cannot be
changed anymore after the first GPU-using command has been issued, for
example ``lbfluid``. If you do not choose the GPU manually before that,
CUDA internally chooses one, which is normally the most powerful GPU
available, but load-independent.::
    
    system=espressomd.System()

    dev=system.cu()
    system.cu(dev)

The first invocation in the sample above return the id of the set graphics card, the second one sets the 
device id.

.. _Creating bonds when particles collide:

Creating bonds when particles collide
-------------------------------------


Please cite :cite:`espresso2` when using dynamic bonding.

With the help of this feature, bonds between particles can be created
automatically during the simulation, every time two particles collide.
This is useful for simulations of chemical reactions and irreversible
adhesion processes. Both, sliding and non-sliding contacts can be created.

The collision detection is controlled via the :attr:`espressomd.system.System.collision_detection` attribute, which is an instance of the class :class:`espressomd.collision_detection.CollisionDetection`.

Several modes are available for different types of binding.

* "bind_centers": adds a pair-bond between two particles at their first collision. By making the bonded interaction `stiff` enough, the particles can be held together after the collision. Note that the particles can still slide on each others' surface, as the pair bond is not directional. This mode is set up as follows::
    import espressomd
    from espressomd.interactions import HarmonicBond
    
    system=espressomd.System()
    bond_centers=HarmonicBond(k=1000,r_0=<CUTOFF>)
    system.bonded_inter.add(bond_centers)
    system.collision_detection.set_params(mode="bind_centers",distance=<CUTOFF>, bond_centers=bond_centers)
  
  The parameters are as follows:
  
    * `distance` is the distance between two particles at which the binding is triggered. This cutoff distance, `<CUTOFF>` in the example above, is typically chosen slightly larger than the particle diameter. It is also a good choice for the equilibrium length of the bond.
    * `bond_centers` is the bonded interaction (an instance of :class:espressomd.interactions.BondedInteraction`) to be created between the particles. No guarantees are made regarding which of the two colliding particles gets the bond. Once there is a bond of this type on any of the colliding particles, no further binding occurs for this pair of particles.

* "bind_at_point_of_collision": this mode prevents sliding of the colliding particles at the contact. This is achieved by
  creating two virtual sites at the point of collision. They are
  rigidly connected to the colliding particles, respectively. A bond is
  then created between the virtual sites, or an angular bond between
  the two colliding particles and the virtual particles. In the latter case,
  the virtual particles are the centers of the angle potentials
  (particle 2 in the description of the angle potential (see :ref:`Bond-angle interactions`).
  Due to the rigid connection between each of the
  particles in the collision and its respective virtual site, a sliding
  at the contact point is no longer possible. See the documentation on
  :ref:`Rigid arrangements of particles` for details. In addition to the bond between the virtual
  sites, the bond between the colliding particles is also created, i.e., the "bind_at_point_of_collision" mode implicitly includes the "bind_centers" mode. You
  can either use a real bonded interaction to prevent wobbling around
  the point of contact or you can use :class:`espressomd.interactions.Virtual` which acts as a marker, only.
  The method is setup as follows::
     
     system.collision_detection.set_params(mode="bind_at_point_of_collision", distance=<CUTOFF>, bond_centers=<BOND_CENTERS>, bond_vs=<BOND_VS>, part_type_vs=<PART_TYPE_VS>, vs_placement=<VS_PLACEMENT>)

  
  The parameters `distance` and `bond_centers` have the same meaning as in the `bind_centers` mode. The remaining parameters are as follows:
    
    * `bond_vs` is the bond to be added between the two virtual sites created on collision. This is either a pair-bond with an equilibrium length matching the distance between the virtual sites, or an angle bond fully stretched in its equilibrium configuration.
    * `part_type_vs` is the particle type assigned to the virtual sites created on collision. In nearly all cases, no non-bonded interactions should be defined for this particle type.
    * `vs_placement` controls, where on the line connecting the centers of the colliding particles, the virtual sites are placed. A value of 0 means that the virtual sites are placed at the same position as the colliding particles on which they are based. A value of 0.5 will result in the virtual sites being placed ad the mid-point between the two colliding particles. A value of 1 will result the virtual site associated to the first colliding particle to be placed at the position of the second colliding particle. In most cases, 0.5, is a good choice. Then, the bond connecting the virtual sites should have an equilibrium length of zero.

* "glue_to_surface": This mode is used to irreversibly attach small particles to the surface of a big particle. It is asymmetric in that several small particles can be bound to a big particle but not vice versa. The small particles can change type after collision to make them `inert`. On collision, a single virtual site is placed and related to the big particle. Then, a bond (`bond_centers`) connects the big and the small particle. A second bond (`bond_vs`) connects the virtual site and the small particle. Further required parameters are:
  
  * `part_type_to_attach_vs_to`: Type of the particle to which the virtual site is attached, i.e., the `big` particle.
  * `part_type_to_be_glued`: Type of the particle bound to the virtual site (the `small` particle).
  * `part_type_after_glueing`: The type assigned to the particle bound to the virtual site (`small` particle) after the collision.
  * `part_type_vs`: Particle type assigned to the virtual site created during the collision. 
  * `distance_glued_particle_to_vs`: Distance of the virtual site to the particle being bound to it (`small` particle).




- "bind_three_particles" allows for the creation of agglomerates which maintain their shape
  similarly to those create by the mode "bind_at_point_of_collision". The present approach works
  without virtual sites. Instead, for each two-particle collision, the
  surrounding is searched for a third particle. If one is found,
  angular bonds are placed to maintain the local shape.
  If all three particles are within the cutoff distance, an angle bond is added
  on each of the three particles in addition
  to the distance based bonds between the particle centers. 
  If two particles are within the cutoff of a centrla particle (e.g., chain of three particles)
  an angle bond is placed on the central particle.
  The angular bonds being added are determined from the angle between the particles.
  This method does not depend on the particles’ rotational
  degrees of freedom being integrated. Virtual sites are also not
  required.
  The method, along with the corresponding bonds are setup as follows::
        
        n_anlge_bonds=181 # 0 to 180 degrees in one degree steps
        for i in range(0,res,1):
           self.s.bonded_inter[i]=Angle_Harmonic(bend=1,phi0=float(i)/(res-1)*np.pi)
        
        # Create the bond passed to bond_centers here and add it to the system
        
        self.s.collision_detection.set_params(mode="bind_three_particles",bond_centers=<BOND_CENTERS>,bond_three_particles=0,three_particle_binding_angle_resolution=res,distance=<CUTOFF>)

  Important: The bonds for the angles are mapped via their numerical bond ids. In this example, ids from 0 to 180 are used. All other bonds required for the simulation need to be added to the system after those bonds. In particular, this applies to the bonded interaction passed via `bond_centers` 


The following limitations currently apply for the collision detection:
* No distinction is currently made between different particle types for the `bind_centers` method.
* The “bind at point of collision” and "glue to surface"  approaches require the feature `VIRTUAL_SITES_RELATIVE` to be activated in `myconfig.hpp`.

* The “bind at point of collision” approach cannot handle collisions
  between virtual sites

Catalytic Reactions
-------------------

With the help of the feature ``CATALYTIC_REACTIONS``, one can define three particle types to act as reactant (e.g. :math:`\mathrm{H_2 O_2}`), catalyzer (e.g. platinum), and product (e.g. :math:`\mathrm{O_2}` and :math:`\mathrm{H_2 O}`). The current setup allows one to simulate active swimmers and their chemical propulsion.

For a Janus swimmer consisting of platinum on one hemisphere and gold on the other hemisphere, both surfaces catalytically induce a reaction. We assume an initial abundance of hydrogen peroxide and absence of products, so that back (recombination) reactions seldomly occur at the surface. A typical model for the propulsion of such a particle assumes

.. math::

    \begin{aligned}
      \mathrm{H_2 O_2} &\xrightarrow{\text{Pt}} \mathrm{2 H^{+} + 2 e^{-} + O_2} \\
      \mathrm{2 H^{+} + 2 e^{-} + H_2 O_2} &\xrightarrow{\text{Au}} \mathrm{2 H_2 O}
    \end{aligned}

That is, catalytic surfaces induce a reactions that produce charged species by consuming hydrogen peroxide. It is the change in distribution of charged species that leads to motion of the swimmer, a process referred to as self-electrophoresis. A minimal model for this would be

.. math::

    \begin{aligned}
      A &\xrightarrow{C^{+}} B \\
      B &\xrightarrow{C^{-}} A
    \end{aligned}

where on the upper half of the catalyst :math:`C^{+}` a species :math:`A` is converted into :math:`B`, and on the lower half :math:`C^{-}` the opposite reaction takes place. Note that when :math:`A` and :math:`B` are charged, this reaction conserves charge, provided the rates are equal.

In |es| the orientation of a catalyzer particle is used to define hemispheres; half spaces going through the particle's center. The reaction region is bounded by the *reaction range*: :math:`r`. Inside the reaction range, we react only rectant-product pairs. The particles in a pair are swapped from hemisphere to another with a rate prescribed by

.. math::

    P_{\text{move}} = 1 - \mathrm{e}^{-k_{\mathrm{ct}}\,\Delta t} ,

with the reaction rate :math:`k_{\mathrm{ct}}` and the simulation time step :math:`\Delta t`. A pair may be swapped only once per MD time step, to avoid a no-net-effect situation. That is, we allow an exchange move only when the following conditions are met:

1. Both partners of the reactant-product pair have to reside within the reaction range.
2. The product has to reside in the upper half-space of the reaction range.
3. The reactant has to reside in the lower half-space of the reaction range.

Self-propulsion is achieved by imposing an interaction asymmetry between the partners of a swapped pair. That is, the heterogeneous distribution of chemical species induced by the swapping leads to a net force on the particle, counter balanced by friction.

To set up the system for catalytic reactions the class :class:`espressomd.reaction.Reaction`
can be used.::

    from espressomd.reaction import Reaction

    system = espressomd.System()

    # setting up particles etc

    r = Reaction(product_type=1, reactant_type=2, catalyzer_type=0, ct_range=2, ct_rate=0.2, eq_rate=0)
    r.start()
    r.stop()

    print r

* the first invocation of ``Reaction``, in the above example,  defines a
  reaction with particles of type number 2 as reactant, type 0 as catalyzer and
  type 1 as product [#1]_. The catalytic reaction rate constant is given by :math:`\mathrm{ct\_rate}`
  [#2]_ and to override the default rate constant for the equilibrium reaction
  ( = 0), one can specify it by as ``eq_rata``.  By default each reactant particle is checked
  against each catalyst particle (``react_once=False``). However, when creating
  smooth surfaces using many catalyst particles, it can be desirable to let the
  reaction rate be independent of the surface density of these particles. That
  is, each particle has a likelihood of reacting in the vicinity of the surface
  (distance is less than :math:`r`) as specified by the rate constant, i.e.,
  *not* according to :math:`P_{\text{cvt}} = 1 - \exp \left( - n k\Delta t
  \right)`, with :math:`n` the number of local catalysts. To accomplish this,
  each reactant is considered only once each time step by using the option
  ``react_once=True`` . The reaction command is set up such that the different
  properties may be influenced individually.

*  ``r.stop()`` disables the reaction. Note that at the moment, there can
   only be one reaction in the simulation.

*  ``print r``  returns the current reaction parameters.

In future versions of |es| the capabilities of the ``CATALYTIC_REACTIONS`` feature may be generalized
to handle multiple reactant, catalyzer, and product types, as well as
more general reaction schemes. Other changes may involve merging the
current implementation with the ``COLLISION_DETECTION`` feature.

.. _galilei_transform: 

Galilei Transform and Particle Velocity Manipulation
----------------------------------------------------

The following class :class:`espressomd.galilei.GalileiTransform` may be useful
in effecting the velocity of the system.::
    
    system = espressomd.System()
    gt = system.galilei

Particle motion and rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    gt.kill_particle_motion()

This command halts all particles in the current simulation, setting
their velocities to zero, as well as their angular momentum if the
option ``rotation`` is specified and the feature ROTATION has been
compiled in.

Forces and torques acting on the particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    gt.kill_particle_forces()

This command sets all forces on the particles to zero, as well as all
torques if the option ``torque`` is specified and the feature ROTATION
has been compiled in.

The centre of mass of the system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    gt.system_CMS()

Returns the center of mass of the whole system. It currently does not
factor in the density fluctuations of the Lattice-Boltzmann fluid.

The centre-of-mass velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
    
    gt.system_CMS_velocity()

Returns the velocity of the center of mass of the whole system.

The Galilei transform
~~~~~~~~~~~~~~~~~~~~~

::

    gt.galilei_transform()

Subtracts the velocity of the center of mass of the whole system from
every particle’s velocity, thereby performing a Galilei transform into
the reference frame of the center of mass of the system. This
transformation is useful for example in combination with the DPD
thermostat, since there, a drift in the velocity of the whole system
leads to an offset in the reported temperature.

.. rubric:: Footnotes

.. [#1]
   Only one type of particle can be assigned to each of these three
   reaction species and no particle type may be assigned to multiple
   species. That is, currently does not support particles of type 1 and
   2 both to be reactants, nor can particles of type 1 be a reactant as
   well as a catalyst. Moreover, only one of these reactions can be
   implemented in a single Tcl script. If, for instance, there is a
   reaction involving particle types 1, 2, and 4, there cannot be a
   second reaction involving particles of type 5, 6, and 8. It is
   however possible to modify the reaction properties for a given set of
   types during the simulation.

.. [#2]
   Currently only strictly positive values of the catalytic conversion
   rate constant are allowed. Setting the value to zero is equivalent to
   ``r.stop()``.
