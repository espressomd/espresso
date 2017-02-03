Running the simulation
======================

``integrate``: Running the simulation
-------------------------------------

To run the integrator use call the
:meth:`espressomd.integrate.Integrator.run` method of the
:class:`espressomd._system.integrator` instance::

    system.integrator.run(steps=number_of_steps)

where ``number_of_steps`` is the number of time steps the integrator
should perform.

|es| uses the Velocity Verlet algorithm for the integration of the equations
of motion.

Note that this implementation of the Velocity Verlet algorithm reuses
forces, that is, they are computed once in the middle of the time step,
but used twice, at the beginning and end. However, in the first time
step after setting up, there are no forces present yet. Therefore, has
to compute them before the first time step. That has two consequences:
first, random forces are redrawn, resulting in a narrower distribution
of the random forces, which we compensate by stretching. Second,
coupling forces of e.g. the Lattice Boltzmann fluid cannot be computed
and are therefore lacking in the first half time step. In order to
minimize these effects, has a quite conservative heuristics to decide
whether a change makes it necessary to recompute forces before the first
time step. Therefore, calling hundred times
:meth:`espressomd.integrate.Integrator.run` with ``steps=1`` does the
same as with ``steps=100``, apart from some small calling overhead.

However, for checkpointing, there is no way for to tell that the forces
that you read back in actually match the parameters that are set.
Therefore, would recompute the forces before the first time step, which
makes it essentially impossible to checkpoint LB simulations, where it
is vital to keep the coupling forces. To work around this, there is
an additional parameter, which tells integrate to not recalculate
the forces for the first time step, but use that the values still stored
with the particles. Use this only if you are absolutely sure that the
forces stored match your current setup!

The opposite problem occurs when timing interactions: In this case, one
would like to recompute the forces, despite the fact that they are
already correctly calculated. To this aim, the option can be used to
enforce force recalculation.

``minimize_energy``: Run steepest descent minimization
------------------------------------------------------

In Python the ``minimize_energy`` functionality can be imported from
:mod:`espressomd.minimize_energy` as class
:class:`espressomd.minimize_energy.MinimizeEnergy`. Alternatively it
is already part of the :class:`espressomd._system.System` class object
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

while the maximal force is bigger than or for at most times. The energy
is relaxed by , while the change per coordinate per step is limited to .
The combination of and can be used to get an poor man’s adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the ``fix`` command or the
``ROTATION_PER_PARTICLE`` feature are not altered.

``change_volume_and_rescale_particles``: Changing the box volume
----------------------------------------------------------------

This is implemented in
:meth:`espressomd._system.System.change_volume_and_rescale_particles`
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
:class:`espressomd._system.System.galilei` has two functions:

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
either the large :attr:`espressomd._system.System.time_step` or the
other :attr:`espressomd._system.System.smaller_time_step`. Setting
:attr:`espressomd._system.System.smaller_time_step` to a positive
value turns on the multi-timestepping algorithm. The ratio
:attr:`espressomd._system.System.time_step`/:attr:`espressomd._system.System.smaller_time_step`
*must* be an integer. Beads are by default associated with
:attr:`espressomd._system.System.time_step`, corresponding to the
particle property
:attr:`espressomd.particle_data.ParticleHandle.smaller_timestep` set
to 0. Setting to
:attr:`espressomd.particle_data.ParticleHandle.smaller_timestep` to 1
associates the particle to the
:attr:`espressomd._system.System.smaller_time_step` integration. The
integrator can be used in the NVE ensemble, as well as with the
Langevin thermostat and the modified Andersen barostat for NVT and NPT
simulations, respectively. See :cite:`bereau15` for more details.
