.. _Running the simulation:

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
coupling forces of e.g. the lattice Boltzmann fluid cannot be computed
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

.. _Run steepest descent minimization:

Run steepest descent minimization
---------------------------------

:func:`espressomd.integrate.Integrator.set_steepest_descent`



This feature is used to propagate each particle by a small distance parallel to the force acting on it.
When only conservative forces for which a potential exists are in use, this is equivalent to a steepest descent energy minimization.
A common application is removing overlap between randomly placed particles.

Please note that the behavior is undefined if a thermostat is activated.
It runs a simple steepest descent algorithm:

Iterate

.. math:: p_i = p_i + \min(\texttt{gamma} \times F_i, \texttt{max_displacement}),

while the maximal force is bigger than ``f_max`` or for at most ``max_steps`` times. The energy
is relaxed by ``gamma``, while the change per coordinate per step is limited to ``max_displacement``.
The combination of ``gamma`` and ``max_displacement`` can be used to get a poor man's adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the ``fix`` and ``rotation`` attribute of particles are not altered.

Usage example::

        system.integrator.set_steepest_descent(
            f_max=0, gamma=0.1, max_displacement=0.1)
        system.integrator.run(20)
        system.integrator.set_vv()  # to switch back to velocity verlet




.. _Integrating rotational degrees of freedom:

Integrating rotational degrees of freedom
-----------------------------------------
When the feature ``ROTATION`` is compiled in, Particles not only have a position, but also an orientation.
Just as a force on a particle leads to an increase in linear velocity, a torque on a particle leads to an increase in angular velocity. The rotational degrees of freedom are also integrated using a velocity Verlet scheme.
When the Langevin thermostat is enabled, the rotational degrees of freedom are also thermalized.

Whether or not rotational degrees of freedom are propagated, is controlled on a per-particle and per-axis level, where the axes are the Cartesian axes of the particle in its body-fixed frame.
It is important to note that starting from version 4.0 and unlike in earlier versions of |es|, the particles' rotation is disabled by default.
In this way, just compiling in the ``ROTATION`` feature no longer changes the physics of the system.

The rotation of a particle is controlled via the :attr:`espressomd.particle_data.ParticleHandle.rotation` property. E.g., the following code adds a particle with rotation on the x axis enabled::

    import espressomd
    s = espressomd.System()
    s.part.add(pos=(0, 0, 0), rotation=(1, 0, 0))

Notes:

* The orientation of a particle is stored as a quaternion in the :attr:`espressomd.particle_data.ParticleHandle.quat` property. For a value of (1,0,0,0), the body and space frames coincide.
* The space-frame direction of the particle's z-axis in its body frame is accessible through the ``espressomd.particle_data.ParticleHandle.director`` property.
* Any other vector can be converted from body to space fixed frame using the ``espressomd.particle_data.ParticleHandle.convert_vector_body_to_space`` method.
* When ``DIPOLES`` are compiled in, the particles dipole moment is always co-aligned with the z-axis in the body-fixed frame.
* Changing the particles dipole moment and director will re-orient the particle such that its z-axis in space frame is aligned parallel to the given vector. No guarantees are made for the other two axes after setting the director or the dipole moment.


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

