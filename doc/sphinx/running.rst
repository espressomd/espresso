.. _Running the simulation:

Running the simulation
======================

.. _Particle integration and propagation:

Particle integration and propagation
------------------------------------

The main integration scheme of |es| is the velocity Verlet algorithm.
A steepest descent algorithm is used to minimize the system.

Additional integration schemes are available, which can be coupled to
thermostats to enable Langevin dynamics, Brownian dynamics, Stokesian dynamics,
dissipative particle dynamics, and simulations in the NpT ensemble.

.. _Integrators:

Integrators
-----------

To run the integrator call the method
:meth:`system.integrate.run() <espressomd.integrate.Integrator.run>`::

    system.integrator.run(number_of_steps, recalc_forces=False, reuse_forces=False)

where ``number_of_steps`` is the number of time steps the integrator should perform.

.. _Velocity Verlet Algorithm:

Velocity Verlet algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^

:meth:`espressomd.integrate.IntegratorHandle.set_vv`

The equations of motion for the trajectory of point-like particles read

.. math:: \dot v_i(t) = F_i(\{x_j\},v_i,t)/m_i \\ \dot x_i(t) = v_i(t),

where :math:`x_i`, :math:`v_i`, :math:`m_i` are position, velocity and mass of
particle :math:`i` and :math:`F_i(\{x_j\},v_i,t)` the forces acting on it.
These forces comprise all interactions with other particles and external fields
as well as non-deterministic contributions described in :ref:`Thermostats`.

For numerical integration, this equation is discretized to the following steps (:cite:`rapaport04a` eqs. 3.5.8 - 3.5.10):

1. Calculate the velocity at the half step

   .. math:: v(t+dt/2) = v(t) + \frac{F(x(t),v(t-dt/2),t)}{m} dt/2

2. Calculate the new position

   .. math:: x(t+dt) = x(t) + v(t+dt/2) dt

3. Calculate the force based on the new position

   .. math:: F = F(x(t+dt), v(t+dt/2), t+dt)

4. Calculate the new velocity

   .. math:: v(t+dt) = v(t+dt/2) + \frac{F(x(t+dt),t+dt)}{m} dt/2

Note that this implementation of the velocity Verlet algorithm reuses
forces in step 1. That is, they are computed once in step 3,
but used twice, in step 4 and in step 1 of the next iteration. In the first time
step after setting up, there are no forces present yet. Therefore, |es| has
to compute them before the first time step. That has two consequences:
first, random forces are redrawn, resulting in a narrower distribution
of the random forces, which we compensate by stretching. Second,
coupling forces of e.g. the lattice-Boltzmann fluid cannot be computed
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

.. _Isotropic NpT integrator:

Isotropic NpT integrator
^^^^^^^^^^^^^^^^^^^^^^^^

:meth:`espressomd.integrate.IntegratorHandle.set_isotropic_npt`

As the NpT thermostat alters the way the equations of motion are integrated, it is
discussed here and only a brief summary is given in :ref:`Thermostats`.

To activate the NpT integrator, use :meth:`~espressomd.integrate.IntegratorHandle.set_isotropic_npt`
with parameters:

* ``ext_pressure``: The external pressure
* ``piston``: The mass of the applied piston
* ``direction``: Flags to enable/disable box dimensions to be subject to fluctuations. By default, all directions are enabled.

Additionally, a NpT thermostat has to be set by :meth:`~espressomd.thermostat.Thermostat.set_npt()`
with parameters:

* ``kT``: Thermal energy of the heat bath
* ``gamma0``: Friction coefficient of the bath
* ``gammav``: Artificial friction coefficient for the volume fluctuations.

A code snippet would look like::

    import espressomd

    system = espressomd.System(box_l=[1, 1, 1])
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0, seed=42)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

The physical meaning of these parameters is described below:

The relaxation towards a desired pressure :math:`P` (parameter ``ext_pressure``)
is enabled by treating the box
volume :math:`V` as a degree of freedom with corresponding momentum :math:`\Pi = Q\dot{V}`,
where :math:`Q` (parameter ``piston``) is an artificial piston mass.
Which box dimensions are affected to change the volume can be controlled by a list of
boolean flags for parameter ``direction``.
An additional energy :math:`H_V = 1/(2Q)\Pi + PV`
associated with the volume is postulated. This results in a "force" on the box such that

.. math:: \dot{\Pi} = \mathcal{P} - P

where

.. math:: \mathcal{P} = \frac{1}{Vd} \sum_{i,j} f_{ij}x_{ij} + \frac{1}{Vd} \sum_i m_i v_i^2

Here :math:`\mathcal{P}` is the instantaneous pressure, :math:`d` the dimension
of the system (number of flags set by ``direction``), :math:`f_{ij}` the
short range interaction force between particles :math:`i` and :math:`j` and
:math:`x_{ij}= x_j - x_i`.

In addition to this deterministic force, a friction :math:`-\frac{\gamma^V}{Q}\Pi(t)`
and noise :math:`\sqrt{k_B T \gamma^V} \eta(t)` are added for the box
volume dynamics and the particle dynamics. This introduces three new parameters:
The friction coefficient for the box :math:`\gamma^V` (parameter ``gammav``),
the friction coefficient of the particles :math:`\gamma^0` (parameter ``gamma0``)
and the thermal energy :math:`k_BT` (parameter ``kT``).
For a discussion of these terms and their discretisation, see :ref:`Langevin thermostat`,
which uses the same approach, but only for particles.
As a result of box geometry changes, the particle positions and velocities have to be rescaled
during integration.

The discretisation consists of the following steps (see :cite:`kolb99a` for a full derivation of the algorithm):

1. Calculate the particle velocities at the half step

   .. math:: v'(t+dt/2) = v(t) + \frac{F(x(t),v(t-dt/2),t)}{m} dt/2

2. Calculate the instantaneous pressure and "volume momentum"

   .. math:: \mathcal{P} = \mathcal{P}(x(t),V(t),f(x(t)), v'(t+dt/2))
   .. math:: \Pi(t+dt/2) = \Pi(t) + (\mathcal{P}-P) dt/2 -\frac{\gamma^V}{Q}\Pi(t) dt/2  +  \sqrt{k_B T \gamma^V dt} \overline{\eta}

3. Calculate box volume and scaling parameter :math:`L` at half step and full step, scale the simulation box accordingly

   .. math:: V(t+dt/2) = V(t) + \frac{\Pi(t+dt/2)}{Q} dt/2
   .. math:: L(t+dt/2) = V(t+dt/2)^{1/d}
   .. math:: V(t+dt) = V(t+dt/2) + \frac{\Pi(t+dt/2)}{Q} dt/2
   .. math:: L(t+dt) = V(t+dt)^{1/d}

4. Update particle positions and scale velocities

   .. math:: x(t+dt) = \frac{L(t+dt)}{L(t)} \left[ x(t) + \frac{L^2(t)}{L^2(t+dt/2)} v(t+dt/2) dt \right]
   .. math:: v(t+dt/2) = \frac{L(t)}{L(t+dt)} v'(t+dt/2)

5. Calculate forces, instantaneous pressure and "volume momentum"

   .. math:: F = F(x(t+dt),v(t+dt/2),t)
   .. math:: \mathcal{P} = \mathcal{P}(x(t+dt),V(t+dt),f(x(t+dt)), v(t+dt/2))
   .. math:: \Pi(t+dt) = \Pi(t+dt/2) + (\mathcal{P}-P) dt/2 -\frac{\gamma^V}{Q}\Pi(t+dt/2) dt/2  +  \sqrt{k_B T \gamma^V dt} \overline{\eta}

   with uncorrelated numbers :math:`\overline{\eta}` drawn from a random uniform process :math:`\eta(t)`

6. Update the velocities

   .. math:: v(t+dt) = v(t+dt/2) + \frac{F(t+dt)}{m} dt/2

Notes:

* The NpT algorithm is only tested for all 3 directions enabled for scaling. Usage of ``direction`` is considered an experimental feature.
* In step 4, only those coordinates are scaled for which ``direction`` is set.
* For the instantaneous pressure, the same limitations of applicability hold as described in :ref:`Pressure`.
* The particle forces :math:`F` include interactions as well as a friction (:math:`\gamma^0`) and noise term (:math:`\sqrt{k_B T \gamma^0 dt} \overline{\eta}`) analogous to the terms in the :ref:`Langevin thermostat`.
* The particle forces are only calculated in step 5 and then reused in step 1 of the next iteration. See :ref:`Velocity Verlet Algorithm` for the implications of that.

.. _Steepest descent:

Steepest descent
^^^^^^^^^^^^^^^^

:meth:`espressomd.integrate.IntegratorHandle.set_steepest_descent`

This feature is used to propagate each particle by a small distance parallel to the force acting on it.
When only conservative forces for which a potential exists are in use, this is equivalent to a steepest descent energy minimization.
A common application is removing overlap between randomly placed particles.

Please note that the behavior is undefined if a thermostat is activated,
in which case the integrator will generate an error. The integrator runs
the following steepest descent algorithm:

.. math:: \vec{r}_{i+1} = \vec{r}_i + \min(\gamma \vec{F}_i, \vec{r}_{\text{max_displacement}}),

while the maximal force/torque is bigger than ``f_max`` or for at most ``steps`` times. The energy
is relaxed by ``gamma``, while the change per coordinate per step is limited to ``max_displacement``.
The combination of ``gamma`` and ``max_displacement`` can be used to get a poor man's adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle,
with ``max_displacement`` interpreted as the maximal rotation angle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the ``fix`` and ``rotation`` attribute of particles are not altered.

Usage example::

    system.integrator.set_steepest_descent(
        f_max=0, gamma=0.1, max_displacement=0.1)
    system.integrator.run(20)   # maximal number of steps
    system.integrator.set_vv()  # to switch back to velocity Verlet

Using a custom convergence criterion
""""""""""""""""""""""""""""""""""""

The ``f_max`` parameter can be set to zero to prevent the integrator from
halting when a specific force/torque is reached. The integration can then
be carried out in a loop with a custom convergence criterion::

    min_sigma = 1  # size of the smallest particle
    max_sigma = 5  # size of the largest particle
    min_dist = 0.0
    system.integrator.set_steepest_descent(f_max=0, gamma=10,
                                           max_displacement=min_sigma * 0.01)
    # gradient descent until particles are separated by at least max_sigma
    while min_dist < max_sigma:
        min_dist = system.analysis.min_dist()
        system.integrator.run(10)
    system.integrator.set_vv()

When writing a custom convergence criterion based on forces or torques, keep
in mind that particles whose motion and rotation are fixed in space along
some or all axes with ``fix`` or ``rotation`` need to be filtered from the
force/torque observable used in the custom convergence criterion. Since these
two properties can be cast to boolean values, they can be used as masks to
remove forces/torques that are ignored by the integrator::

    particles = system.part.all()
    max_force = np.max(np.linalg.norm(particles.f * np.logical_not(particles.fix), axis=1))
    max_torque = np.max(np.linalg.norm(particles.torque_lab * np.logical_not(particles.rotation), axis=1))

Virtual sites can also be an issue since the force on the virtual site is
transferred to the target particle at the beginning of the integration loop.
The correct forces need to be re-calculated after running the integration::

    def convergence_criterion(forces):
        '''Function that decides when the gradient descent has converged'''
        return ...
    p1 = system.part.add(pos=[0, 0, 0], type=1)
    p2 = system.part.add(pos=[0, 0, 0.1], type=1)
    p2.vs_auto_relate_to(p1)
    system.integrator.set_steepest_descent(f_max=800, gamma=1.0, max_displacement=0.01)
    while convergence_criterion(system.part.all().f):
        system.integrator.run(10)
        system.integrator.run(0, recalc_forces=True)  # re-calculate forces from virtual sites
    system.integrator.set_vv()

The algorithm can also be used for energy minimization::

    # minimize until energy difference < 5% or energy < 1e-3
    system.integrator.set_steepest_descent(f_max=0, gamma=1.0, max_displacement=0.01)
    relative_energy_change = float('inf')
    relative_energy_change_threshold = 0.05
    energy_threshold = 1e-3
    energy_old = system.analysis.energy()['total']
    print(f'Energy: {energy_old:.2e}')
    for i in range(20):
        system.integrator.run(50)
        energy = system.analysis.energy()['total']
        print(f'Energy: {energy:.2e}')
        relative_energy_change = (energy_old - energy) / energy_old
        if relative_energy_change < relative_energy_change_threshold or energy < energy_threshold:
            break
        energy_old = energy
    else:
        print(f'Energy minimization did not converge in {i + 1} cycles')
    system.integrator.set_vv()

Please note that not all features support energy calculation.

.. _Brownian Dynamics:

Brownian Dynamics
^^^^^^^^^^^^^^^^^

Brownian Dynamics integrator :cite:`schlick10a`.
See details in :ref:`Brownian thermostat`.

.. _Stokesian Dynamics:

Stokesian Dynamics
^^^^^^^^^^^^^^^^^^

.. note::

    Requires ``STOKESIAN_DYNAMICS`` external feature, enabled with
    ``-DWITH_STOKESIAN_DYNAMICS=ON``.

:meth:`espressomd.integrate.IntegratorHandle.set_stokesian_dynamics`

The Stokesian Dynamics method allows to study the behaviour of spherical
particles in a viscous fluid. It is targeted at systems with very low Reynolds
numbers. In such systems, particles come to a rest almost immediately as soon as
any force on them is removed. In other words, motion has no memory of the past.

The integration scheme is relatively simple. Only the particles' positions,
radii and forces (including torques) are needed to compute the momentary
velocities (including angular velocities). The particle positions are
integrated by the simple Euler scheme.

The computation of the velocities is an approximation with good results
in the far field.
The Stokesian Dynamics method is only available for open systems,
i.e. no periodic boundary conditions are supported. The box size has
no effect either.

The Stokesian Dynamics method is outlined in :cite:`durlofsky87a`.

The following minimal example illustrates how to use the SDM in |es|::

    import espressomd
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.periodicity = [False, False, False]
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.part.add(pos=[0, 0, 0], rotation=[1, 0, 0])
    system.integrator.set_stokesian_dynamics(viscosity=1.0, radii={0: 1.0})
    system.integrator.run(100)

Because there is no force on the particle yet, nothing will move. You will need
to add your own actors to the system. The parameter ``radii`` is a dictionary
that maps particle types to different radii. ``viscosity`` is the dynamic
viscosity of the ambient infinite fluid. There are additional optional
parameters for ``set_stokesian_dynamics()``. For more information, see
:py:meth:`espressomd.integrate.IntegratorHandle.set_stokesian_dynamics()`.

Note that this setup represents a system at zero temperature. In order to
thermalize the system, the SD thermostat needs to be activated (see
:ref:`Stokesian thermostat`).

.. _Important_SD:

Important
"""""""""

The particles must be prevented from overlapping. It is mathematically allowed
for the particles to overlap to a certain degree. However, once the distance
of the sphere centers is less than 2/3 of the sphere diameter, the mobility
matrix is no longer positive definite and the Stokesian Dynamics integrator
will fail. Therefore, the particle centers must be kept apart from each
other by a strongly repulsive potential, for example the WCA potential
that is set to the appropriate particle radius (for more information about
the available interaction types see :ref:`Non-bonded interactions`).

The current implementation of SD only includes the far field approximation.
The near field (so-called lubrication) correction is planned. For now,
Stokesian Dynamics provides a good approximation of the hydrodynamics
in dilute systems where the average distance between particles is several
sphere diameters.


.. _Thermostats:

Thermostats
-----------

To add a thermostat, call the appropriate setter::

    system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=41)

The different thermostats available in |es| will be described in the following
subsections.

You may combine different thermostats at your own risk by turning them on
one by one. The list of active thermostats can be cleared at any time with
:py:meth:`system.thermostat.turn_off() <espressomd.thermostat.Thermostat.turn_off>`.
Not all combinations of thermostats are allowed, though (see
:py:func:`espressomd.thermostat.AssertThermostatType` for details).
Some integrators only work with a specific thermostat and throw an
error otherwise. Note that there is only one temperature for all
thermostats, although for some thermostats like the Langevin thermostat,
particles can be assigned individual temperatures.

Since |es| does not enforce a particular unit system, it cannot know about
the current value of the Boltzmann constant. Therefore, when specifying
the temperature of a thermostat, you actually do not define the
temperature, but the value of the thermal energy :math:`k_B T` in the
current unit system (see the discussion on units, Section :ref:`On units`).

All thermostats have a ``seed`` argument that controls the state of the random
number generator (Philox Counter-based RNG). This seed is required on first
activation of a thermostat, unless stated otherwise. It can be omitted in
subsequent calls of the method that activates the same thermostat. The random
sequence also depends on the thermostats counters that are
incremented after each integration step.

.. _Langevin thermostat:

Langevin thermostat
^^^^^^^^^^^^^^^^^^^

In order to activate the Langevin thermostat the member function
:py:meth:`~espressomd.thermostat.Thermostat.set_langevin` of the thermostat
class :class:`espressomd.thermostat.Thermostat` has to be invoked.
Best explained in an example::

    import espressomd
    system = espressomd.System(box_l=[1, 1, 1])
    system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=41)

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
This reduces the accuracy of the velocity Verlet integrator
by one order in :math:`dt` because forces are now velocity-dependent.

The random process :math:`\eta(t)` is discretized by drawing an uncorrelated random number
:math:`\overline{\eta}` for each component of all the particle forces.
The distribution of :math:`\overline{\eta}` is uniform and satisfies

.. math:: <\overline{\eta}> = 0 , <\overline{\eta}\overline{\eta}> = 1/dt

If the feature ``ROTATION`` is compiled in, the rotational degrees of freedom are
also coupled to the thermostat. If only the first two arguments are
specified then the friction coefficient for the rotation is set to the
same value as that for the translation.
A separate rotational friction coefficient can be set by inputting
``gamma_rotate``. The two options allow one to switch the translational and rotational
thermalization on or off separately, maintaining the frictional behavior. This
can be useful, for instance, in high Péclet number active matter systems, where
one wants to thermalize only the rotational degrees of freedom while
translational degrees of freedom are affected by the self-propulsion.

The keywords ``gamma`` and ``gamma_rotate`` can be specified as a scalar,
or, with feature ``PARTICLE_ANISOTROPY`` compiled in, as the three eigenvalues
of the respective friction coefficient tensor. This is enables the simulation of
the anisotropic diffusion of anisotropic colloids (rods, etc.).

Using the Langevin thermostat, it is possible to set a temperature and a
friction coefficient for every particle individually via the feature
``THERMOSTAT_PER_PARTICLE``.  Consult the reference of the ``part`` command
(chapter :ref:`Setting up particles`) for information on how to achieve this.

.. _Brownian thermostat:

Brownian thermostat
^^^^^^^^^^^^^^^^^^^

Brownian thermostat is a formal name of a thermostat enabling the
Brownian Dynamics feature (see :cite:`schlick10a`) which implies
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
    system = espressomd.System(box_l=[1, 1, 1])
    system.thermostat.set_brownian(kT=1.0, gamma=1.0, seed=41)
    system.integrator.set_brownian_dynamics()

where ``gamma`` (hereinafter :math:`\gamma`) is a viscous friction coefficient.
In terms of the Python interface and setup, the Brownian thermostat is very
similar to the :ref:`Langevin thermostat`. The feature
``THERMOSTAT_PER_PARTICLE`` is used to control the per-particle
temperature and the friction coefficient setup. The major differences are
its internal integrator implementation and other temporal constraints.
The integrator is still a symplectic velocity Verlet-like one.
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

Note: the velocity random walk is propagated from zero at each step.

A rotational motion is implemented similarly.
Note: the rotational Brownian dynamics implementation is compatible with particles which have
the isotropic moment of inertia tensor only. Otherwise, the viscous terminal angular velocity
is not defined, i.e. it has no constant direction over the time.

.. _Isotropic NpT thermostat:

Isotropic NpT thermostat
^^^^^^^^^^^^^^^^^^^^^^^^

This feature allows to simulate an (on average) homogeneous and isotropic system in the NpT ensemble.
In order to use this feature, ``NPT`` has to be defined in the :file:`myconfig.hpp`.
Activate the NpT thermostat with the command :py:meth:`~espressomd.thermostat.Thermostat.set_npt`
and setup the integrator for the NpT ensemble with :py:meth:`~espressomd.integrate.IntegratorHandle.set_isotropic_npt`.

For example::

    import espressomd

    system = espressomd.System(box_l=[1, 1, 1])
    system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=1.0, seed=41)
    system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

For an explanation of the algorithm involved, see :ref:`Isotropic NpT integrator`.

Be aware that this feature is neither properly examined for all systems
nor is it maintained regularly. If you use it and notice strange
behavior, please contribute to solving the problem.

.. _Dissipative Particle Dynamics (DPD):

Dissipative Particle Dynamics (DPD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The DPD thermostat adds friction and noise to the particle
dynamics like the :ref:`Langevin thermostat`, but these
are not applied to every particle individually but instead
encoded in a dissipative interaction between particles :cite:`soddeman03a`.

To realize a complete DPD fluid model in |es|, three parts are needed:
the DPD thermostat, which controls the temperate, a dissipative interaction
between the particles that make up the fluid, see :ref:`DPD interaction`,
and a repulsive conservative force, see :ref:`Hat interaction`.

The temperature is set via
:py:meth:`espressomd.thermostat.Thermostat.set_dpd`
which takes ``kT`` and ``seed`` as arguments.

The friction coefficients and cutoff are controlled via the
:ref:`DPD interaction` on a per type-pair basis.

The friction (dissipative) and noise (random) term are coupled via the
fluctuation-dissipation theorem. The friction term is a function of the
relative velocity of particle pairs. The DPD thermostat is better for
dynamics than the Langevin thermostat, since it mimics hydrodynamics in
the system.

As a conservative force any interaction potential can be used,
see :ref:`Isotropic non-bonded interactions`. A common choice is
a force ramp which is implemented as :ref:`Hat interaction`.

A complete example of setting up a DPD fluid and running it
to sample the equation of state can be found in :file:`/samples/dpd.py`.

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

.. _LB thermostat:

Lattice-Boltzmann thermostat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`Lattice-Boltzmann` thermostat acts similar to the :ref:`Langevin thermostat` in that the governing equation for particles is

.. math::  m_i \dot{v}_i(t) = f_i(\{x_j\},v_i,t) - \gamma (v_i(t)-u(x_i(t),t)) + \sqrt{2\gamma k_B T} \eta_i(t).

where :math:`u(x,t)` is the fluid velocity at position :math:`x` and time :math:`t`.
To preserve momentum, an equal and opposite friction force and random force act on the fluid.

Numerically the fluid velocity is determined from the lattice-Boltzmann node velocities
by interpolating as described in :ref:`Interpolating velocities`.
The backcoupling of friction forces and noise to the fluid is also done by distributing those forces amongst the nearest LB nodes.
Details for both the interpolation and the force distribution can be found in :cite:`ahlrichs99a` and :cite:`duenweg08a`.

The LB fluid can be used to thermalize particles, while also including their hydrodynamic interactions.
The LB thermostat expects an instance of either :class:`espressomd.lb.LBFluid` or :class:`espressomd.lb.LBFluidGPU`.
Temperature is set via the ``kT`` argument of the LB fluid.

The magnitude of the frictional coupling can be adjusted by the
parameter ``gamma``. To enable the LB thermostat, use::

    import espressomd
    import espressomd.lb
    system = espressomd.System(box_l=[1, 1, 1])
    lbf = espressomd.lb.LBFluid(agrid=1, dens=1, visc=1, tau=0.01)
    system.actors.add(lbf)
    system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

No other thermostatting mechanism is necessary
then. Please switch off any other thermostat before starting the LB
thermostatting mechanism.

The LBM implementation provides a fully thermalized LB fluid, all
nonconserved modes, including the pressure tensor, fluctuate correctly
according to the given temperature and the relaxation parameters. All
fluctuations can be switched off by setting the temperature to 0.

.. note:: Coupling between LB and MD only happens if the LB thermostat is set with a :math:`\gamma \ge 0.0`.

.. _Stokesian thermostat:

Stokesian thermostat
^^^^^^^^^^^^^^^^^^^^

.. note::

    Requires ``STOKESIAN_DYNAMICS`` external feature, enabled with
    ``-DWITH_STOKESIAN_DYNAMICS=ON``.

In order to thermalize a Stokesian Dynamics simulation, the SD thermostat
needs to be activated via::

    import espressomd
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.periodicity = [False, False, False]
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.part.add(pos=[0, 0, 0], rotation=[1, 0, 0], ext_force=[0, 0, -1])
    system.thermostat.set_stokesian(kT=1.0, seed=43)
    system.integrator.set_stokesian_dynamics(viscosity=1.0, radii={0: 1.0})
    system.integrator.run(100)

where ``kT`` denotes the desired temperature of the system, and ``seed`` the
seed for the random number generator.
