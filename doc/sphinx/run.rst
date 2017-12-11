Running the simulation
======================

.. _Integrator:

Integrator
----------

To run the integrator call the method
:meth:`espressomd.integrate.Integrator.run`::

    system.integrator.run(steps=number_of_steps, recalc_forces=False, reuse_forces=False)

where ``number_of_steps`` is the number of time steps the integrator
should perform.

|es| uses the Velocity Verlet algorithm for the integration of the equations
of motion.

Note that this implementation of the Velocity Verlet algorithm reuses
forces. That is, they are computed once in the middle of the time step,
but used twice, at the beginning and end. In the first time
step after setting up, there are no forces present yet. Therefore, |es| has
to compute them before the first time step. That has two consequences:
first, random forces are redrawn, resulting in a narrower distribution
of the random forces, which we compensate by stretching. Second,
coupling forces of e.g. the Lattice Boltzmann fluid cannot be computed
and are therefore lacking in the first half time step. In order to
minimize these effects, |es| has a quite conservative heuristics to decide
whether a change makes it necessary to recompute forces before the first
time step. Therefore, calling 100 times
:meth:`espressomd.integrate.Integrator.run` with ``steps=1`` does the
same as with ``steps=100``, apart from some small calling overhead.

However, for checkpointing, there is no way for |es| to tell that the forces
that you read back in actually match the parameters that are set.
Therefore, |es| would recompute the forces before the first time step, which
makes it essentially impossible to checkpoint LB simulations, where it
is vital to keep the coupling forces. To work around this, there is
an additional parameter ``reuse_forces``, which tells integrate to not recalculate
the forces for the first time step, but use that the values still stored
with the particles. Use this only if you are absolutely sure that the
forces stored match your current setup!

The opposite problem occurs when timing interactions: In this case, one
would like to recompute the forces, despite the fact that they are
already correctly calculated. To this aim, the option ``recalc_forces`` can be used to
enforce force recalculation.

Run steepest descent minimization
---------------------------------

In Python the ``minimize_energy`` functionality can be imported from
:mod:`espressomd.minimize_energy` as class
:class:`espressomd.minimize_energy.MinimizeEnergy`. Alternatively it
is already part of the :class:`espressomd.system.System` class object
and can be called from there (second variant)::

    espressomd.minimize_energy.init(
        f_max = <double>,
        gamma = <double>,
        max_steps = <double>,
        max_displacement = <double>)
    espressomd.minimize_energy.minimize()

    system.minimize_energy.init(
        f_max = <double>,
        gamma = <double>,
        max_steps = <double>,
        max_displacement = <double>)
    system.minimize_energy.minimize()

This command runs a steepest descent energy minimization on the system.
Please note that the behaviour is undefined if either a thermostat,
Maggs electrostatics or Lattice-Boltzmann is activated. It runs a simple
steepest descent algorithm:

Iterate

.. math:: p_i = p_i + \min(\texttt{gamma} \times F_i, \texttt{max_displacement}),

while the maximal force is bigger than ``f_max`` or for at most ``max_steps`` times. The energy
is relaxed by ``gamma``, while the change per coordinate per step is limited to ``max_displacement``.
The combination of ``gamma`` and ``max_displacement`` can be used to get a poor man’s adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the ``fix`` command or the
``ROTATION_PER_PARTICLE`` feature are not altered.

Changing the box volume
-----------------------

This is implemented in
:meth:`espressomd.system.System.change_volume_and_rescale_particles`
with the parameters ``d_new`` for the new length and ``dir`` for the
coordinates to work on and ``"xyz"`` for isotropic.

Changes the volume of either a cubic simulation box to the new volume or
its given x-/y-/z-/xyz-extension to the new box-length, and
isotropically adjusts the particles coordinates as well. The function
returns the new volume of the deformed simulation box.

Stopping particles
------------------

To stop particles you can use the functionality implemented in the
:mod:`espressomd.galilei` module.  The corresponding class
:class:`espressomd.galilei.GalileiTransform` which is wrapped inside
the :class:`espressomd.system.System` instance as
:class:`espressomd.system.System.galilei` has two functions:

- :meth:`espressomd.galilei.GalileiTransform.kill_particle_motion`:
   halts all particles in the current simulation, setting their
   velocities to zero, as well as their angular momentum if the
   feature ``ROTATION`` has been compiled in.

- :meth:`espressomd.galilei.GalileiTransform.kill_particle_forces`:
   sets all forces on the particles to zero, as well as all torques if
   the feature ``ROTATION`` has been compiled in.

Multi-timestepping
------------------

Required feature: ``MULTI_TIMESTEP``

The multi-timestepping integrator allows to run two concurrent
integration time steps within a simulation, associating beads with
either the large :attr:`espressomd.system.System.time_step` or the
other :attr:`espressomd.system.System.smaller_time_step`. Setting
:attr:`espressomd.system.System.smaller_time_step` to a positive
value turns on the multi-timestepping algorithm. Beads are by default associated with
:attr:`espressomd.system.System.time_step`, corresponding to the
particle property
:attr:`espressomd.particle_data.ParticleHandle.smaller_timestep` set
to 0. Setting to
:attr:`espressomd.particle_data.ParticleHandle.smaller_timestep` to 1
associates the particle to the
:attr:`espressomd.system.System.smaller_time_step` integration. The
integrator can be used in the NVE ensemble, as well as with the
Langevin thermostat and the modified Andersen barostat for NVT and NPT
simulations, respectively. See :cite:`bereau15` for more details.

Reaction Ensemble
-----------------

.. note:: The whole Reaction Ensemble module uses Monte Carlo moves which require potential energies. Therefore the Reaction Ensemble requires support for energy calculations for all interactions which are used in the simulation.

The reaction ensemble :cite:`smith94a,turner2008simulation` allows to simulate
chemical reactions which can be represented by the general equation:

.. math::

   \mathrm{\nu_1 S_1 +\ \dots\  \nu_l S_l\ \rightleftharpoons\ \nu_m S_m +\ \dots\ \nu_z S_z }
       \label{general-eq}

where :math:`\nu_i` is the stoichiometric coefficient of species
:math:`S_i`. By convention, stoichiometric coefficents of the
species on the left-hand side of the reaction (*reactants*) attain
negative values, and those on the right-hand side (*products*) attain
positive values, so that the reaction can be equivalently written as

.. math::

   \mathrm{\sum_i \nu_i S_i = 0} \,.
       \label{general-eq-sum}


The equilibrium constant of the reaction is then given as

.. math::

   K = \exp(-\Delta_{\mathrm{r}}G^{\ominus} / k_B T)
       \quad\text{with}\quad
       \Delta_{\mathrm{r}}G^{\ominus} = \sum_i \nu_i \mu_i^{\ominus}\,.
       \label{Keq}


Here :math:`k_B` is the Boltzmann constant, :math:`T` is temperature,
:math:`\Delta_{\mathrm{r}}G^{\ominus}` standard Gibbs free energy change
of the reaction, and :math:`\mu_i^{\ominus}` the standard chemical
potential (per particle) of species :math:`i`. Note that thermodynamic equilibrium is
independent of the direction in which we write the reaction. If it is
written with left and righ-hand side swapped, 
both :math:`\Delta_{\mathrm{r}}G^{\ominus}` and the stoichiometric
coefficients attain opposite signs, and the equilibrium constant attains the inverse value. 
Further, note that the equilibrium constant :math:`K` is the
dimensionless *thermodynamic, concentration-based* equilibrium constant,
defined as

.. math::

   K(c^{\ominus}) = (c^{\ominus})^{-\bar\nu} \prod_i (c_i)^{\nu_i}

wher :math:`\bar\nu=\sum_i \nu_i`, and :math:`c^{\ominus}` is the reference concentration,
at which the standard chemical potential :math:`\Delta_{\mathrm{r}}G^{\ominus}` was determined.
In practice, this constant is often used with the dimension of :math:`(c^{\ominus})^{\bar\nu}`

.. math::

   K_c(c^{\ominus}) = K(c^{\ominus})\times (c^{\ominus})^{\bar\nu}

A simulation in
the reaction ensemble consists of two types of moves: the *reaction move*
and the *configuration move*. The configuration move changes the configuration
of the system. It is not performed by the Reaction Ensemble module, and can be
performed by a suitable molecular dynamics or a Monte Carlo scheme. The
``reacton_ensemble`` command takes care only of the reaction moves.
In the *forward* reaction, the appropriate number of reactants (given by
:math:`\nu_i`) is removed from the system, and the concomitant number of
products is inserted into the system. In the *backward* reaction,
reactants and products exchange their roles. The acceptance probability
:math:`P^{\xi}` for move from state :math:`o` to :math:`n` reaction
ensemble is given by the criterion :cite:`smith94a`

.. math::

   P^{\xi} = \text{min}\biggl(1,V^{\bar\nu\xi}\Gamma^{\xi}e^{-\beta\Delta E}\prod_{i=1}\frac{N_i^0!}{(N_i^0+\nu_{i}\xi)!}
       \label{eq:Pacc}
       \biggr),

where :math:`\Delta E=E_\mathrm{new}-E_\mathrm{old}` is the change in potential energy,
:math:`V` is the simulation box volume,
and :math:`\beta=1/k_\mathrm{B}T`. 
The extent of reaction, :math:`\xi=1` for the forward, and
:math:`\xi=-1` for the backward direction. 
The parameter :math:`\Gamma` proportional to the reaction constant. It is defined as

.. math::

   \Gamma = \prod_i \Bigl(\frac{\left<N_i\right>}{V} \Bigr)^{\bar\nu} = V^{-\bar\nu} \prod_i \left<N_i\right>^{\nu_i} = K_c(c^{\ominus}=1/\sigma^3)

where :math:`\left<N_i\right>/V` is the average number density of particles of type :math:`i`.
Note that the dimension of :math:`\Gamma` is :math:`V^{\bar\nu}`, therefore its
units must be consistent with the units in which Espresso measures the box volume,
i.e. :math:`\sigma^3`.
   
It is often convenient, and in some cases even necessary, that some particles
representing reactants are not removed from or placed at randomly in the system
but their identity is changed to that of the products, or vice versa in the
backward direction.  A typical example is the ionization reaction of weak
polyelectrolytes, where the ionizable groups on the polymer have to remain on
the polymer chain after the reaction.  The replacement rule is that the identity of a given reactant type is
changed to the corresponding product type as long as the corresponding
coefficients allow for it.  Corresponding means having the same position (index) in
the python lists of reactants and products which are used to set up the
reaction.
For a description of the available methods see :mod:`espressomd.reaction_ensemble`

Converting tabulated reaction constants to internal units in Espresso 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The implementation in Espresso requires that the dimension of :math:`\Gamma` 
is consistent with the internal unit of volume, :math:`\sigma^3`.
The tabulated values of equilibrium constants for reactions in solution, :math:`K_c`, typically use
:math:`c^{\ominus} = 1\,\mathrm{moldm^{-3}}` as the reference concentration, 
and have the dimension of :math:`(c^{\ominus})^{\bar\nu}`.  To be used with Espresso, the
value of :math:`K_c` has to be converted as

.. math::

   \Gamma = K_c(c^{\ominus} = 1/\sigma^3) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}}) 
   \Bigl( N_{\mathrm{A}}\bigl(\frac{\sigma}{\mathrm{dm}}\bigr)^3\Bigr)^{\bar\nu}
   
where :math:`N_{\mathrm{A}}` is the Avogardo number.  For gas-phase reactions,
the pressure-based eaction constant, :math:`K_p` is often used, which can
be converted to :math:`K_c` as

.. math::

   K_p(p^{\ominus}=1\,\mathrm{atm}) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}}) \biggl(\frac{c^{\ominus}RT}{p^{\ominus}}\biggr)^{\bar\nu},

where :math:`p^{\ominus}=1\,\mathrm{atm}` is the standard pressure. 



.. The text below is commented-out because it is still an open research question how it should be used correctly.
..
.. This can be used to include water autoprotolysis in the implicit solvent simulation, 
.. by means of a reaction:
.. 
.. .. math::
.. 
..    \mathrm{2 H_2O \rightleftharpoons\ H_3O^+ + OH^- } \,,
.. 
.. 
.. add the following ex nihilo reactions to Espresso. (:math:`\emptyset`, read ex
.. nihilo). Ex nihilo means that the reaction has no reactants or products.
.. Therefore, if :math:`\emptyset` is a product, particles vanish and if
.. :math:`\emptyset` is an reactant, then particles are created ex nihilo:
.. 
.. .. math::
.. 
..    \mathrm{\emptyset \rightleftharpoons\ H_3O^+ + OH^- }  \,, 
.. 
.. with reaction constant K
.. 
.. .. math::
.. 
..    \mathrm{H_3O^+ + OH^- \rightleftharpoons\ \emptyset} \,, 
.. 
.. with reaction constant 1/K. K is given implicitly as a function of the apparent dissociation
.. constant :math:`K_w=10^{-14} \rm{mol^2/l^2}=x\cdot \rm{1/(\sigma^3)^2}` such that the dimensionless is
.. :math:`K=(x\cdot \rm{1/(\sigma^3)^2})/(\beta P^0)^{\overline{\nu}}` with
.. :math:`\overline{\nu}=2` for the dissociation reaction and where x is
.. the value of the apparent dissociation constant that is converted from
.. :math:`\rm{mol^2/l^2}` to a number density in :math:`1/(\sigma^3)^2`,
.. where :math:`\sigma` is the simulation length unit. If :math:`\beta` and
.. :math:`P^0` are provided in simulation units this will make :math:`K`
.. dimensionless. As a test for the autodissociation of water a big
.. simulation box can be set up and the autodissociation reaction can be
.. performed. Then the box should fill with the correct number of protons
.. and hydroxide ions (check for the number of protons and hydroxide ions
.. in the given simulation volume and compare this to the expected value at
.. pH 7). Further the :math:`pK_w=14` should be reproduced -also in the
.. case of an initial excess of acid or base in the simulation box. Note
.. that this only works for big enough volumes.




Wang-Landau Reaction Ensemble 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. .. note:: Requires support for energy calculations for all used interactions since it uses Monte-Carlo moves which use energies in one way or the other.

Combination of the Reaction Ensemble with the Wang-Landau algorithm 
:cite:`wang01a`
allows for enhanced sampling of the reacting system, and
and for the determination of the density of states with respect 
to the reaction coordinate or with respect to some other collective
variable :cite:`landsgesell16a`. Here the 1/t Wang-Landau
algorithm :cite:`belardinelli07a` is implemented since it
does not suffer from systematic errors. Additionally to the above
commands for the reaction ensemble use the following commands for the
Wang-Landau reaction ensemble. For a description of the available methods see :mod:`espressomd.reaction_ensemble`:

Constant pH simulation using the Reaction Ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. .. note:: Requires support for energy calculations for all used interactions since it uses Monte-Carlo moves which use energies.

In the constant pH method due to Reed and Reed
:cite:`reed92a` it is possible to set the chemical potential
of :math:`H^{+}` ions, assuming that the simulated system is coupled to an
infinite reservoir. This value is the used to simulate dissociation
equilibrium of acids and bases. Under certain conditions, the constant
pH method can yield equivalent results as the reaction ensemble :cite:`landsgesell16b`. However, it
treats the chemical potential of :math:`H^{+}` ions and their actual
number in the simulation box as independent variables, which can lead to
serious artifacts. 
The constant pH method can be used within the reaction ensemble module by
initializing the reactions with the standard commands of the reaction ensemble. 

The dissociation constant, which is the input of the constant pH method, is the equilibrium
constant :math:`K_c` for the following reaction:

.. math::

   \mathrm{HA \rightleftharpoons\ H^+ + A^- } \,,

For an example of how to setup
a Constant pH simulation, see the file in the testsuite directory. 
For a description of the available methods see :mod:`espressomd.reaction_ensemble`:

Grand canonical ensemble simulation using the Reaction Ensemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a special case, all stoichiometric coefficients on one side of the chemical
reaction can be set to zero.  Such reaction creates particles *ex nihilo*, and 
is equivalent to exchange with a reservoir. Then the simulation in the reaction ensemble becomes equivalent with the
grandcanonical simulation. Formally, this can be expressed by the reaction

.. math::
 
    \mathrm{\emptyset \rightleftharpoons\ \nu_A A  }  \,, 

where, if :math:`\nu_A=1`, the reaction constant :math:`\Gamma` defines the chemical potential of species A.
However, if :math:`\nu_A\neq 1`, the statistics of the reaction ensemble becomes
equivalent to the grandcanonical only in the limit of large average number of species A in the box.
It the reaction contains more than one product, then the reaction constant
:math:`\Gamma` defines only the sum of their chemical potentials but not the
chemical potential of each product alone.

.. Since the Reaction Ensemble acceptance transition probability can be
.. derived from the grand canonical acceptance transition probability we
.. can use the reaction ensemble to implement grand canonical simulation
.. moves. This is done via adding reactions that only have reactants (for the
.. deletion of particles) or only have products (for the creation of
.. particles). There exists a one to one mapping of the expressions in the
.. grand canonical transition probabilities and the expressions in the
.. reaction ensemble transition probabilities.



Integrating rotational degrees of freedom
-----------------------------------------
When the feature ROTATION is compiled in, Particles not only have a position, but also an orientation.
Just as a force on a particle leads to an increase in linear velocity, a torque on a particle leads to an increase in angular velocity. The rotational degrees of freedom are also integrated using a velocity Verlet scheme.
When the Langevin thermostat is enabled, the rotational degrees of freedom are also thermalized.

Whether or not rotational degrees of freedom are propagated, is controlled on a per-particle and per-axis level, where the axes are the Cartesian axes of the particle in its body-fixed frame.
It is important to note that starting from version 4.0 and unlike in earlier versions of |es|, the particles' rotation is disabled by default.
In this way, just compiling in the ROTATION feature no longer changes the physics of the system.

The rotation of a particle is controlled via the :attr:`espressomd.particle_data.ParticleHandle.rotation` property. E.g., the following code adds a particle with rotation on the x axis enabled:::
    
    import espressomd
    s=espressomd.System()
    s.part.add(pos=(0,0,0),rotation=(1,0,0))

Notes:

* The orientation of a particle is stored as a quaternion in the :attr:`espressomd.particle_data.ParticleHandle.quat` property. For a value of (1,0,0,0), the body and space frames coincide. 
* The space-frame direction of the particle's z-axis in its body frame is accessible through the `espressomd.particle_data.ParticleHandle.director` property.
* Any other vector can be converted from body to space fixed frame using the `espressomd.particle_data.ParticleHandle.convert_vector_body_to_space` method.
* When DIPOLES are compiled in, the particles dipole moment is always co-aligned with the z-axis in the body-fixed frame.
* Changing the particles dipole moment and director will re-orient the particle such that its z-axis in space frame is aligned parallel to the given vector. No guarantees are made for the other two axes after setting the direcotr or the dipole moment.


The following particle properties are related to rotation:

* :attr:`espressomd.particle_data.ParticleHandle.dip`
* :attr:`espressomd.particle_data.ParticleHandle.director`
* :attr:`espressomd.particle_data.ParticleHandle.ext_torque`
* :attr:`espressomd.particle_data.ParticleHandle.gamma_rot`
* :attr:`espressomd.particle_data.ParticleHandle.gamma_rot`
* :attr:`espressomd.particle_data.ParticleHandle.omega_body`
* :attr:`espressomd.particle_data.ParticleHandle.omega_lab`
* :attr:`espressomd.particle_data.ParticleHandle.quat`
* :attr:`espressomd.particle_data.ParticleHandle.rinertia`
* :attr:`espressomd.particle_data.ParticleHandle.rotation`
* :attr:`espressomd.particle_data.ParticleHandle.torque_lab`





