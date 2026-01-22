.. _Magnetostatics:

Magnetostatics
==============

.. _Dipolar interaction:

Dipolar interaction
-------------------

|es| contains methods to calculate the interactions between point dipoles

.. math::

   U^{Dip}(\vec{r}) = D \cdot \left( \frac{(\vec{\mu}_i \cdot \vec{\mu}_j)}{r^3}
     - \frac{3  (\vec{\mu}_i \cdot \vec{r})  (\vec{\mu}_j \cdot \vec{r}) }{r^5} \right)

where :math:`r=|\vec{r}|`.
The prefactor :math:`D` is can be set by the user and is given by

.. math::
   D =\frac{\mu_0 \mu}{4\pi}
   :label: dipolar_prefactor

where :math:`\mu_0` and :math:`\mu` are the vacuum permittivity and the
relative permittivity of the background material, respectively.

Magnetostatic interactions are activated when attached to the system::

    import espressomd
    import espressomd.magnetostatics

    system = espressomd.System(box_l=[10, 10, 10])
    system.time_step = 0.01
    system.part.add(pos=[[0, 0, 0], [1, 1, 1]], dip=2 * [(1, 0, 0)],
                    rotation=2 * [(True, True, True)])

    actor = espressomd.magnetostatics.DipolarDirectSum(prefactor=1.)
    system.magnetostatics.solver = actor

The solver can be detached with either::

    system.magnetostatics.solver = None

or::

    system.magnetostatics.clear()


.. _Dipolar P3M:

Dipolar P3M
~~~~~~~~~~~

:class:`espressomd.magnetostatics.DipolarP3M`

.. note::

    Requires feature ``DIPOLES`` and
    external feature ``FFTW``, enabled with
    ``-D ESPRESSO_BUILD_WITH_FFTW=ON``.

This is the dipolar version of the P3M algorithm, described in :cite:`cerda08d`.

Make sure that you know the relevance of the P3M parameters before using
P3M! If you are not sure, read the following references:
:cite:`ewald21a,hockney88a,kolafa92a,deserno98a,deserno98b,deserno00e,deserno00b,cerda08d`.

Note that dipolar P3M does not work with non-cubic boxes.


The parameters of the dipolar P3M method can be tuned automatically, by
providing ``accuracy=<TARGET_ACCURACY>`` to the method. It is also possible to
pass a subset of the method parameters such as ``mesh``. In that case, only
the omitted parameters are tuned::

    import espressomd.magnetostatics as magnetostatics
    p3m = magnetostatics.DipolarP3M(prefactor=1, mesh=32, accuracy=1E-4)
    system.magnetostatics.solver = p3m

It is important to note that the error estimates given in :cite:`cerda08d`
used in the tuning contain assumptions about the system. In particular, a
homogeneous system is assumed. If this is no longer the case during the
simulation, actual force and torque errors can be significantly larger.


.. _Dipolar Layer Correction (DLC):

Dipolar Layer Correction (DLC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`espressomd.magnetostatics.DLC`

.. note::

    Requires feature ``DIPOLES``.

The dipolar layer correction (DLC) is used in conjunction with the dipolar P3M
method to calculate dipolar interactions in a 2D-periodic system.
It is based on :cite:`brodka04a` and the dipolar version of
:ref:`Electrostatic Layer Correction (ELC)`.

Usage notes:

* The non-periodic direction is always the **z-direction**.

* The method relies on a slab of the simulation box perpendicular to the
  z-direction not to contain particles. The size in z-direction of this slab
  is controlled by the ``gap_size`` parameter. The user has to ensure that
  no particles enter this region by means of constraints or by fixing the
  particles' z-coordinate. When particles enter the slab of the specified
  size, an error will be thrown.

* The method can be tuned using the ``accuracy`` parameter. In contrast to
  the electrostatic method, it refers to the energy. Furthermore, it is
  assumed that all dipole moment are as large as the largest of the dipoles
  in the system.

* When the base solver is not a P3M method, metallic epsilon is assumed.

The method is used as follows::

    import espressomd.magnetostatics
    dp3m = espressomd.magnetostatics.DipolarP3M(prefactor=1, accuracy=1E-4)
    mdlc = espressomd.magnetostatics.DLC(actor=dp3m, maxPWerror=1E-5, gap_size=2.)
    system.magnetostatics.solver = mdlc


.. _Dipolar direct sum:

Dipolar direct sum
------------------

:class:`~espressomd.magnetostatics.DipolarDirectSum`

.. note::

    Requires feature ``DIPOLES`` and optionally
    external feature ``CUDA``, enabled with
    ``-D ESPRESSO_BUILD_WITH_CUDA=ON``.

This interaction calculates energies and forces between dipoles by
explicitly summing over all pairs. For the directions in which the
system is periodic (as defined by ``system.periodicity``), it applies the
minimum image convention, i.e. the interaction is effectively cut off at
half a box length. This can be changed by adding replicas.

The direct summation method is mainly intended for non-periodic systems
which cannot be solved using the dipolar P3M method.
Due to the long-range nature of dipolar interactions, direct summation with
minimum image convention does not yield good accuracy with periodic systems.

The CPU implementation performs the calculation in double-precision,
while the GPU implementation performs the calculations in single-precision.

The only required parameter is the prefactor :eq:`dipolar_prefactor`::

    import espressomd.magnetostatics
    dds = espressomd.magnetostatics.DipolarDirectSum(prefactor=1., n_replicas=2, gpu=False)
    system.magnetostatics.solver = dds

The optional argument ``n_replicas`` adds periodic copies to the system along
periodic directions. In that case, the minimum image convention is no longer used.
Additionally, enabling the ``DIPOLE_FIELD_TRACKING`` feature enables either
implementation to calculate the total dipole field at the position of each
magnetic particle in the primary simulation box. These values are stored in
the particle handle's ``dip_fld`` property and can be accessed directly or
via an observable.

Both the CPU and GPU implementations support MPI-parallelization.
The GPU implementation doesn't support multi-GPU acceleration.


.. _ScaFaCoS magnetostatics:

ScaFaCoS magnetostatics
-----------------------

:class:`espressomd.magnetostatics.Scafacos`

.. note::

    Requires feature ``SCAFACOS_DIPOLES`` and
    external feature ``SCAFACOS``, enabled with
    ``-D ESPRESSO_BUILD_WITH_SCAFACOS=ON``.

|es| can use the methods from the ScaFaCoS *Scalable fast Coulomb solvers*
library for dipoles, if the methods support dipolar calculations.
Dipolar calculations are only included in the ``dipoles`` branch of
the ScaFaCoS code. The specific methods available can be queried with
:meth:`espressomd.electrostatics.Scafacos.get_available_methods`.

To use ScaFaCoS, create an instance of :class:`~espressomd.magnetostatics.Scafacos`
and attach it to the system. Three parameters have to be specified:
``prefactor``, ``method_name``, ``method_params``. The method-specific
parameters are described in the ScaFaCoS manual. In addition, methods
supporting tuning have a parameter ``tolerance_field`` which sets the desired
root mean square accuracy for the magnetic field.

Here is an example with the P2NFFT method for 2D-periodic systems::

    import espressomd.magnetostatics
    scafacos = espressomd.magnetostatics.Scafacos(
        prefactor=1., method_name="p2nfft",
        method_params={
            "p2nfft_verbose_tuning": 0,
            "pnfft_N": "80,80,160",
            "pnfft_window_name": "bspline",
            "pnfft_m": "4",
            "p2nfft_ignore_tolerance": "1",
            "pnfft_diff_ik": "0",
            "p2nfft_r_cut": "6",
            "p2nfft_alpha": "0.8",
            "p2nfft_epsB": "0.05"})
    system.magnetostatics.solver = scafacos

For details of the various methods and their parameters please refer to
the ScaFaCoS manual. To use this feature, ScaFaCoS has to be built as a
shared library.
