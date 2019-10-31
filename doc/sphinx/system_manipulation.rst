.. _System manipulation:

System manipulation
===================

.. _Changing the box volume:

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

.. _Stopping particles:

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

.. _Fixing the center of mass:

Fixing the center of mass
-------------------------

This interaction type applies a constraint on particles of the specified
types such that during the integration the center of mass of these particles is
fixed. This is accomplished as follows: The sum of all the forces acting
on particles of type are calculated. These include all the forces due to
other interaction types and also the thermostat. Next a force equal in
magnitude, but in the opposite direction is applied to all the
particles. This force is divided on the particles of type relative to
their respective mass. Under periodic boundary conditions, this fixes
the itinerant center of mass, that is, the one obtained from the
unfolded coordinates.

Center of mass fixing can be activated via  :class:`espressomd.system.System`::

    system.comfixed.types = list_of_types_to_fix

.. _Capping the force during warmup:

Capping the force during warmup
-------------------------------

Non-bonded interactions are often used to model the hard core repulsion
between particles. Most of the potentials in the section are therefore
singular at zero distance, and forces usually become very large for
distances below the particle size. This is not a problem during the
simulation, as particles will simply avoid overlapping. However,
creating an initial dense random configuration without overlap is often
difficult. By artificially capping the forces, it is possible to simulate a system
with overlaps. By gradually raising the cap value, possible overlaps
become unfavorable, and the system equilibrates to an overlap-free
configuration.

Force capping can be activated via  :class:`espressomd.system.System`::

    system.force_cap = F_max

This command will limit the magnitude of the force to :math:`r F_\mathrm{max}`.
Energies are not affected by the capping, so the energy can be used to
identify the remaining overlap. The force capping is switched off by setting
:math:`F_\mathrm{max}=0`.

.. _Galilei Transform and Particle Velocity Manipulation:

Galilei Transform and Particle Velocity Manipulation
----------------------------------------------------

The following class :class:`espressomd.galilei.GalileiTransform` may be useful
in affecting the velocity of the system. ::

    system = espressomd.System()
    gt = system.galilei

* Particle motion and rotation

::

    gt.kill_particle_motion()

This command halts all particles in the current simulation, setting
their velocities to zero, as well as their angular momentum if the
option ``rotation`` is specified and the feature ``ROTATION`` has been
compiled in.

* Forces and torques acting on the particles

  ::

    gt.kill_particle_forces()

  This command sets all forces on the particles to zero, as well as all
  torques if the option ``torque`` is specified and the feature ``ROTATION``
  has been compiled in.

* The center of mass of the system

  ::

    gt.system_CMS()

  Returns the center of mass of the whole system. It currently does not
  factor in the density fluctuations of the lattice-Boltzmann fluid.

* The center-of-mass velocity

  ::

    gt.system_CMS_velocity()

  Returns the velocity of the center of mass of the whole system.

* The Galilei transform

  ::

    gt.galilei_transform()

  Subtracts the velocity of the center of mass of the whole system from
  every particle's velocity, thereby performing a Galilei transform into
  the reference frame of the center of mass of the system. This
  transformation is useful for example in combination with the DPD
  thermostat, since there, a drift in the velocity of the whole system
  leads to an offset in the reported temperature.

