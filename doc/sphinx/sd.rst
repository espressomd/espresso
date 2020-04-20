.. _Stokesian_Dynamics:

Stokesian Dynamics
==================

The Stokesian Dynamics method allows to study the behaviour of spherical
particles in a viscous fluid. It is targeted at systems with very low Reynolds
numbers. In such systems, particles stop moving immediately as soon as any
force on them is removed. In other words, motion has no memory of the past. 

The integration scheme is relatively simple. Only the particle's positions,
radii and forces (including torques) are needed to compute the momentary
velocities (including angular velocities). The particle positions are
integrated by simple Euler scheme.

..note ::The computation of the velocities is an approximation with good results in the far field and in the very near field.
..note ::The Stokesian Dynamics method is only available for open systems, i.e. no periodic boundary conditions are supported. The box size has no effect either.

.. _Setting up an SD system:

Setting up an SD system:
------------------------

The following minimal example illustrates how to use the SDM in |es|::

    import espressomd
    system = espressomd.System(box_l=[10, 10, 10])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
    system.thermostat.set_sd(viscosity=1.0,
                             device="cpu",
                             radii={0: 1.0},
                             flags=flags.SELF_MOBILITY | flags.PAIR_MOBILITY | flags.FTS)
    system.integrator.set_sd()
    system.integrator.run(100)

Because there is no force on the particle yet, nothing will move. You will need
to add your own actors to the system.
