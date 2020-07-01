.. _Stokesian Dynamics:

Stokesian Dynamics
==================

The Stokesian Dynamics method allows to study the behaviour of spherical
particles in a viscous fluid. It is targeted at systems with very low Reynolds
numbers. In such systems, particles stop moving immediately as soon as any
force on them is removed. In other words, motion has no memory of the past.

The integration scheme is relatively simple. Only the particle's positions,
radii and forces (including torques) are needed to compute the momentary
velocities (including angular velocities). The particle positions are
integrated by the simple Euler scheme.

The Stokesian Dynamics method is outlined in :cite:`durlofsky87a`.

.. note::
  The computation of the velocities is an approximation with good results
  in the far field.

.. note::
  The Stokesian Dynamics method is only available for open systems,
  i.e. no periodic boundary conditions are supported. The box size has
  no effect either.

.. _Setting up an SD system:

Setting up an SD system
-----------------------

The following minimal example illustrates how to use the SDM in |es|::

    import espressomd
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.part.add(pos=[0, 0, 0], rotation=[1, 0, 0])
    system.integrator.set_sd(viscosity=1.0, radii={0: 1.0})

    system.integrator.run(100)

Because there is no force on the particle yet, nothing will move. You will need
to add your own actors to the system. The parameter ``radii`` is a dictionary
that maps particle types to different radii. ``viscosity`` is the dynamic
viscosity of the ambient infinite fluid. There are additional optional
parameters for ``set_sd()``. For more information, see
:py:meth:`espressomd.integrate.IntegratorHandle.set_sd()`.

Note that this setup represents a system at zero temperature. In order to
thermalize the system, the SD thermostat can be activated via an additional
command::

    system.thermostat.set_sd(kT=1.0, seed=43)

``kT`` denotes the desired temperature of the system, and ``seed`` the seed for
the random number generator of the Stokesian Dynamics thermostat. It is
independent from other random number generators inside |es|.

.. _Important_SD:

Important
^^^^^^^^^

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
