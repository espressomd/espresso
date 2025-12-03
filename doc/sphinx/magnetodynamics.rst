.. _Magnetodynamics:

Magnetodynamics
===============

|es| contains methods to simulate the internal magnetisation dynamics
of magnetic particles and the interactions between point dipoles.
The default behaviour in |es| is that, once a dipole moment is assigned
to a particle, it will strictly follow the rotation of the particle’s quaternion.
This is called the **fixed point dipole** model.
This model is predominantly used to simulate magnetic soft matter.

However, this approximation does not always represent realistic behaviour
of a single-domain magnetic nanoparticle — it only holds in the limit of
infinitely high magnetic anisotropy energy.
In reality, several **internal relaxation mechanisms** exist,
so the dipole moment is not necessarily co-aligned with the particle quaternion,
nor does it perfectly follow its motion.

To incorporate this phenomenology in simulations, |es| offers several models, listed below.

.. _Thermal_Stoner_Wohlfarth:

Thermal Stoner–Wohlfarth
------------------------

.. note::

    Requires feature ``THERMAL_STONER_WOHLFARTH`` and
    optionally feature ``DIPOLE_FIELD_TRACKING``, as well as
    external feature ``NLOPT``,
    enabled with ``-D ESPRESSO_BUILD_WITH_NLOPT=ON``.

The **thermal Stoner–Wohlfarth (SW)** model includes Néel relaxation in simulations
of single-domain magnetic nanoparticles.
This is an implementation of the algorithm originally presented in :cite:`mostarac25a`.
If you use this method in your work, please cite the original publication, in addition to |es|.

In interacting systems, the method relies on the ``DIPOLE_FIELD_TRACKING`` feature.
Make sure you use magnetostatics actors that support this feature.

---

The magnetic energy of a Stoner–Wohlfarth particle (at :math:`T=0`) is given by

.. math::

   U = - \mu_0 \mu (\vec{e} \cdot \vec{H}) - K V (\vec{e} \cdot \vec{n})^2

The algorithm can be logically separated into three steps:

**Step 1 — Critical Field and Energy Minima**

First, the **critical field** :math:`h_\mathrm{cr}` is calculated based on
the current angle :math:`\phi` between the magnetic field vector :math:`\vec{h}`
and the particle’s anisotropy axis.

For a given :math:`\phi`, the algorithm finds the angle :math:`\theta'_\mathrm{min}`
that minimizes the total magnetic energy and is closest to the previous dipole-moment state.
This is ensured by initializing the minimizer with the previous particle state.

If the field acting on the particle is less than :math:`h_\mathrm{cr}`,
the algorithm proceeds to find the angles :math:`\theta'_\mathrm{max}` and
:math:`\theta''_\mathrm{max}` that maximize the total magnetic energy on
both sides of :math:`\theta'_\mathrm{min}`.

**Step 2 — Energy Barrier and Néel Relaxation**

Using these extrema, the algorithm calculates the **energy barriers**
on both sides of :math:`\theta'_\mathrm{min}`:

.. math::

   \begin{aligned}
   \Delta E'  &= \frac{1}{K V} \left| U(\phi, \theta'_\mathrm{max})  - U(\phi, \theta'_\mathrm{min}) \right|, \\
   \Delta E'' &= \frac{1}{K V} \left| U(\phi, \theta''_\mathrm{max}) - U(\phi, \theta'_\mathrm{min}) \right|.
   \end{aligned}

The smaller of the two barriers is chosen:

.. math::

   \Delta E = \min(\Delta E', \Delta E'').

This barrier is used to estimate a **characteristic timescale** for the Néel relaxation process:

.. math::

   \tau_N = \frac{\tau_D}{2\sigma} \sqrt{\frac{\pi}{\sigma}} e^{\Delta E \sigma},

and the **transition probability** (neglecting back-switching) is

.. math::

   p = 1 - e^{-\delta t / \tau_N},

where :math:`\delta t` is the integration time step.

**Step 3 — Trial Move and Dipole Update**

A **trial move** is made by drawing a random number and comparing it with
the transition probability, similar to a Metropolis Monte Carlo step.

- If the trial move is **successful**, the algorithm finds a new minimum
  :math:`\theta''_\mathrm{min}` and aligns the dipole moment accordingly.
- Otherwise, the dipole moment remains aligned with :math:`\theta'_\mathrm{min}`.

For more details, particularly if you are unsure about the physical quantities involved,
please refer to the original publication.

The example below shows how to set up and parametrise a particle to be used by
the **thermal Stoner-Wohlfarth** solver. Note that ``is_enabled`` needs to be set
to ``True`` explicitly on the virtual site.
Moreover, the virtual sites that are tagged for magnetodynamics must be set
to use the ``Propagation.ROT_VS_INDEPENDENT`` propagation mode.

.. code-block:: python

    import espressomd
    import espressomd.propagation
    Propagation = espressomd.propagation.Propagation

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.001  # MD time step in simulation units
    system.thermostat.set_langevin(kT=1., gamma=75., gamma_rotation=25., seed=42)

    # One particle with thermal Stoner-Wohlfarth enabled
    p1 = system.part.add(pos=[1, 1, 1])
    p1.director = (1, 0, 0) # easy axis direction
    p1.rotation = (True, True, True)  # allow particle rotation

    p2 = system.part.add(pos=p1.pos)
    p2.dip = (1.75,0,0) # set dipole moment for the virtual particle in reduced units
    p2.rotation = (False, False, False) # disable rotations of the virtual site
    p2.magnetodynamics = {
      'is_enabled': True,
      'anisotropy_field_inv': 0.175, # inverse anisotropy field (1/H_k) in reduced units
      'sat_mag': 1.75, # saturation magnetisation in reduced units
      'anisotropy_energy': 5., # anisotropy energy K * V in reduced units
      'sw_dt_incr': 1.0e-10, # kinetic Monte Carlo time increment [s]
      'sw_tau0_inv': 1.0e9  # inverse attempt time (1/tau_0) [1/s]
      }
    # make virtual and set the proper propagation mode for magnetodynamics
    p2.vs_auto_relate_to(p1)
    p2.propagation = Propagation.TRANS_VS_RELATIVE | Propagation.ROT_VS_INDEPENDENT

