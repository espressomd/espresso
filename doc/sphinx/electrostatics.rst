.. _Electrostatics:

Electrostatics
==============

The Coulomb (or electrostatic) interaction is defined as
follows. For a pair of particles at distance :math:`r` with charges
:math:`q_1` and :math:`q_2`, the interaction is given by

.. math:: U_C(r)=C \cdot \frac{q_1 q_2}{r}

where

.. math::
   C=\frac{1}{4\pi \varepsilon_0 \varepsilon_r}
   :label: coulomb_prefactor

is a prefactor which can be set by the user. The commonly used Bjerrum length
:math:`l_B = e^2 / (4 \pi \varepsilon_0 \varepsilon_r k_B T)` is the length at
which the Coulomb energy between two unit charges is equal to the thermal
energy :math:`k_B T`.
Based on this length, the prefactor is given by :math:`C=l_B k_B T / e^2`.

Computing electrostatic interactions is computationally very expensive.
|es| features some state-of-the-art algorithms to deal with these
interactions as efficiently as possible, but almost all of them require
some knowledge to use them properly. Uneducated use can result in
completely unphysical simulations.

Coulomb interactions have to be added to the list of active actors of the system object to become
active. This is done by calling the add-method of :attr:`espressomd.system.System.actors`.
Only one electrostatics method can be active at any time.

Note that using the electrostatic interaction also requires assigning charges to
the particles via the particle property
:py:attr:`espressomd.particle_data.ParticleHandle.q`.

All solvers need a prefactor and a set of other required parameters.
This example shows the general usage of the electrostatic method ``P3M``.
An instance of the solver is created and added to the actors list, at which
point it will be automatically activated. This activation will internally
call a tuning function to achieve the requested accuracy::

    import espressomd
    import espressomd.electrostatics

    system = espressomd.System(box_l=[10, 10, 10])
    system.time_step = 0.01
    system.part.add(pos=[[0, 0, 0], [1, 1, 1]], q=[-1, 1])
    solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-3)
    system.actors.add(solver)

where the prefactor is defined as :math:`C` in Eqn. :eq:`coulomb_prefactor`.

.. _Coulomb P3M:

Coulomb P3M
-----------

:class:`espressomd.electrostatics.P3M`

For this feature to work, you need to have the ``fftw3`` library
installed on your system. In |es|, you can check if it is compiled in by
checking for the feature ``FFTW`` with ``espressomd.features``.
P3M requires full periodicity (1 1 1). When using a non-metallic dielectric
constant (``epsilon != 0.0``), the box must be cubic.
Make sure that you know the relevance of the P3M parameters before using P3M!
If you are not sure, read the following references:
:cite:`ewald21a,hockney88,kolafa92a,deserno98a,deserno98b,deserno00e,deserno00b,cerda08d`.

.. _Tuning Coulomb P3M:

Tuning Coulomb P3M
~~~~~~~~~~~~~~~~~~

The tuning method is called when the handle of the Coulomb P3M is added to the
actor list. At this point, the system should already contain the charged
particles. Set parameters are fixed and not changed by the tuning algorithm.
This can be useful to speed up the tuning during testing or if the parameters
are already known.

To prevent the automatic tuning, set the ``tune`` parameter to ``False``.
To manually tune or retune P3M, call :meth:`espressomd.electrostatics.P3M.tune
<espressomd.electrostatics.ElectrostaticInteraction.tune>`.
Note, however, that this is a method the P3M object inherited from
:attr:`espressomd.electrostatics.ElectrostaticInteraction`.
All parameters passed to the method are fixed in the tuning routine. If not
specified in the ``tune()`` method, the parameters ``prefactor`` and
``accuracy`` are reused.

It is not easy to calculate the various parameters of the P3M method
such that the method provides the desired accuracy at maximum speed. To
simplify this, it provides a function to automatically tune the algorithm.
Note that for this function to work properly, your system should already
contain an initial configuration of charges and the correct initial box
size. Also note that the provided tuning algorithms works very well on
homogeneous charge distributions, but might not achieve the requested
precision for highly inhomogeneous or symmetric systems. For example,
because of the nature of the P3M algorithm, systems are problematic
where most charges are placed in one plane, one small region, or on a
regular grid.

The function employs the analytical expression of the error estimate for
the P3M method :cite:`hockney88` and its real space error :cite:`kolafa92a` to
obtain sets of parameters that yield the desired accuracy, then it measures how
long it takes to compute the Coulomb interaction using these parameter sets and
chooses the set with the shortest run time.

After execution the tuning routines report the tested parameter sets,
the corresponding k-space and real-space errors and the timings needed
for force calculations. In the output, the timings are given in units of
milliseconds, length scales are in units of inverse box lengths.

.. _Coulomb P3M on GPU:

Coulomb P3M on GPU
~~~~~~~~~~~~~~~~~~

:class:`espressomd.electrostatics.P3MGPU`

The GPU implementation of P3M calculates the far field portion on the GPU.
It uses the same parameters and interface functionality as the CPU version of
the solver. It should be noted that this does not always provide significant
increase in performance. Furthermore it computes the far field interactions
with only single precision which limits the maximum precision. The algorithm
does not work in combination with the electrostatic extensions
:ref:`Dielectric interfaces with the ICC* algorithm <Dielectric interfaces with the ICC algorithm>`
and :ref:`Electrostatic Layer Correction (ELC)`.

.. _Debye-Hückel potential:

Debye-Hückel potential
----------------------

:class:`espressomd.electrostatics.DH`

The Debye-Hückel electrostatic potential is defined by

.. math:: U^{C-DH} = C \cdot \frac{q_1 q_2 \exp(-\kappa r)}{r}\quad \mathrm{for}\quad r<r_{\mathrm{cut}}

where :math:`C` is defined as in Eqn. :eq:`coulomb_prefactor` and
:math:`\kappa` is the inverse Debye screening length.
The Debye-Hückel potential is an approximate method for calculating
electrostatic interactions, but technically it is treated as other
short-ranged non-bonding potentials. For :math:`r > r_{\textrm{cut}}` it is
set to zero which introduces a step in energy. Therefore, it introduces
fluctuations in energy.

For :math:`\kappa = 0`, this corresponds to the plain Coulomb potential.

.. _Reaction Field method:

Reaction Field method
---------------------

:class:`espressomd.electrostatics.ReactionField`

The Reaction Field electrostatic potential is defined by

.. math:: U^{C-RF} = C \cdot q_1 q_2 \left[\frac{1}{r} - \frac{B r^2}{2r_{\mathrm{cut}}^3} - \frac{1 - B/2}{r_{\mathrm{cut}}}\right] \quad \mathrm{for}\quad r<r_{\mathrm{cut}}

where :math:`C` is defined as in Eqn. :eq:`coulomb_prefactor` and :math:`B`
is defined as:

.. math:: B = \frac{2(\varepsilon_1 - \varepsilon_2)(1 + \kappa r_{\mathrm{cut}}) - \varepsilon_2 (\kappa r_{\mathrm{cut}})^2}{(\varepsilon_1 + 2\varepsilon_2)(1 + \kappa r_{\mathrm{cut}}) + \varepsilon_2 (\kappa r_{\mathrm{cut}})^2}

with :math:`\kappa` the inverse Debye screening length, :math:`\varepsilon_1` the dielectric
constant inside the cavity and :math:`\varepsilon_2` the dielectric constant
outside the cavity :cite:`tironi95a`.

The term in :math:`1 - B/2` is a correction to make the
potential continuous at :math:`r = r_{\mathrm{cut}}`.


.. _Dielectric interfaces with the ICC algorithm:

Dielectric interfaces with the ICC\ :math:`\star` algorithm
-----------------------------------------------------------

:class:`espressomd.electrostatic_extensions.ICC`

The ICC\ :math:`\star` algorithm allows to take into account arbitrarily shaped
dielectric interfaces and dynamic charge induction. For instance, it can be
used to simulate a curved metallic boundary. This is done by iterating the
charge on a set of spatially fixed *ICC particles* until they correctly
represent the influence of the dielectric discontinuity. All *ICC particles*
need a certain area, normal vector and dielectric constant to specify the
surface. ICC relies on a Coulomb solver that is already initialized. So far, it
is implemented and well tested with the Coulomb solver P3M. ICC is an |es|
actor and can be activated via::

    import espressomd.electrostatic_extensions
    icc = espressomd.electrostatic_extensions.ICC(...)
    system.actors.add(icc)

The ICC particles are setup as normal |es| particles. Note that they should
be fixed in space and need an initial non-zero charge. The following example
sets up parallel metallic plates and activates ICC::

    # Set the ICC line density and calculate the number of
    # ICC particles according to the box size
    l = 3.2
    nicc = int(box_l / l)
    nicc_per_electrode = nicc * nicc
    nicc_tot = 2 * nicc_per_electrode
    iccArea = box_l * box_l / nicc_per_electrode
    l = box_l / nicc

    # Lists to collect required parameters
    iccNormals = []
    iccAreas = []
    iccSigmas = []
    iccEpsilons = []

    # Add the fixed ICC particles:

    # Left electrode (normal [0,0,1])
    for xi in xrange(nicc):
        for yi in xrange(nicc):
            system.part.add(pos=[l * xi, l * yi, 0], q=-0.0001, fix=3*[True], type=icc_type)
    iccNormals.extend([0, 0, 1] * nicc_per_electrode)

    # Right electrode (normal [0,0,-1])
    for xi in xrange(nicc):
        for yi in xrange(nicc):
            system.part.add(pos=[l * xi, l * yi, box_l], q=0.0001, fix=3*[True], type=icc_type)
    iccNormals.extend([0, 0, -1] * nicc_per_electrode)

    # Common area, sigma and metallic epsilon
    iccAreas.extend([iccArea] * nicc_tot)
    iccSigmas.extend([0] * nicc_tot)
    iccEpsilons.extend([100000] * nicc_tot)

    icc = ICC(first_id=0,
              n_icc=nicc_tot,
              convergence=1e-4,
              relaxation=0.75,
              ext_field=[0, 0, 0],
              max_iterations=100,
              eps_out=1.0,
              normals=iccNormals,
              areas=iccAreas,
              sigmas=iccSigmas,
              epsilons=iccEpsilons)

    system.actors.add(icc)


With each iteration, ICC has to solve electrostatics which can severely slow
down the integration. The performance can be improved by using multiple cores,
a minimal set of ICC particles and convergence and relaxation parameters that
result in a minimal number of iterations. Also please make sure to read the
corresponding articles, mainly :cite:`arnold13a,tyagi10a,kesselheim11a` before
using it.

.. _Electrostatic Layer Correction (ELC):

Electrostatic Layer Correction (ELC)
------------------------------------

:class:`espressomd.electrostatics.ELC`

*ELC* is an extension of the P3M electrostatics solver for explicit 2D periodic
systems. It can account for different dielectric jumps on both sides of the
non-periodic direction. In more detail, it is a special procedure that
converts a 3D electrostatic method to a 2D method in computational order N.
Currently, it only supports P3M without GPU. This means,
that you will first have to set up the P3M algorithm before using ELC.
The periodicity has to be set to (1 1 1). *ELC* cancels the electrostatic
contribution of the periodic replica in **z-direction**. Make sure that you
read the papers on ELC (:cite:`arnold02c,arnold02d,tyagi08a`) before using it.
See :ref:`ELC theory` for more details.

Usage notes:

* The non-periodic direction is always the **z-direction**.

* The method relies on a slab of the simulation box perpendicular to the
  z-direction not to contain particles. The size in z-direction of this slab
  is controlled by the ``gap_size`` parameter. The user has to ensure that
  no particles enter this region by means of constraints or by fixing the
  particles' z-coordinate. When particles enter the slab of the specified
  size, an error will be thrown.

*ELC* is an |es| actor and is used with::

    import espressomd.electrostatics
    p3m = espressomd.electrostatics.P3M(prefactor=1, accuracy=1e-4)
    elc = espressomd.electrostatics.ELC(p3m_actor=p3m, gap_size=box_l * 0.2, maxPWerror=1e-3)
    system.actors.add(elc)

*ELC* can also be used to simulate 2D periodic systems with image charges,
specified by dielectric contrasts on the non-periodic boundaries
(:cite:`tyagi08a`). This is achieved by setting the dielectric jump from the
simulation region (*middle*) to *bottom* (at :math:`z=0`) and from *middle* to
*top* (at :math:`z = L_z - h`), where :math:`L_z` denotes the box length in
:math:`z`-direction and :math:`h` the gap size. The corresponding expressions
are :math:`\Delta_t=\frac{\varepsilon_m-\varepsilon_t}{\varepsilon_m+\varepsilon_t}`
and :math:`\Delta_b=\frac{\varepsilon_m-\varepsilon_b}{\varepsilon_m+\varepsilon_b}`::

    elc = espressomd.electrostatics.ELC(p3m_actor=p3m, gap_size=box_l * 0.2, maxPWerror=1e-3,
                                        delta_mid_top=0.9, delta_mid_bot=0.1)

The fully metallic case :math:`\Delta_t=\Delta_b=-1` would lead to divergence
of the forces/energies in *ELC* and is therefore only possible with the
``const_pot`` option.

Toggle ``const_pot`` on to maintain a constant electric potential difference
``pot_diff`` between the xy-planes at :math:`z=0` and :math:`z = L_z - h`::

    elc = espressomd.electrostatics.ELC(p3m_actor=p3m, gap_size=box_l * 0.2, maxPWerror=1e-3,
                                        const_pot=True, delta_mid_bot=100.0)

This is done by countering the total dipole moment of the system with the
electric field :math:`E_{\textrm{induced}}` and superposing a homogeneous
electric field :math:`E_{\textrm{applied}} = \frac{U}{L}` to retain :math:`U`.
This mimics the induction of surface charges
:math:`\pm\sigma = E_{\textrm{induced}} \cdot \varepsilon_0`
for planar electrodes at :math:`z=0` and :math:`z=L_z - h` in a capacitor
connected to a battery with voltage ``pot_diff``.


.. _MMM1D:

MMM1D
-----

:class:`espressomd.electrostatics.MMM1D`

.. note::
    Required features: ``ELECTROSTATICS`` for MMM1D, the GPU version
    additionally needs the features ``CUDA`` and ``MMM1D_GPU``.

Please cite :cite:`arnold05a` when using MMM1D. See :ref:`MMM1D theory` for
the details.

MMM1D is used with::

    import espressomd.electrostatics
    mmm1d = espressomd.electrostatics.MMM1D(prefactor=C, far_switch_radius=fr,
                                            maxPWerror=err, tune=False, bessel_cutoff=bc)
    mmm1d = espressomd.electrostatics.MMM1D(prefactor=C, maxPWerror=err)

where the prefactor :math:`C` is defined in Eqn. :eq:`coulomb_prefactor`.
MMM1D Coulomb method for systems with periodicity (0 0 1). Needs the
N-squared cell system (see section :ref:`Cellsystems`). The first form sets parameters
manually. The switch radius determines at which xy-distance the force
calculation switches from the near to the far formula. The Bessel cutoff
does not need to be specified as it is automatically determined from the
particle distances and maximal pairwise error. The second tuning form
just takes the maximal pairwise error and tries out a lot of switching
radii to find out the fastest one. If this takes too long, you can
change the value of the ``timings`` argument of the
:class:`~espressomd.electrostatics.MMM1D` class,
which controls the number of test force calculations.

.. _MMM1D on GPU:

MMM1D on GPU
~~~~~~~~~~~~

:class:`espressomd.electrostatics.MMM1DGPU`

MMM1D is also available in a GPU implementation. Unlike its CPU
counterpart, it does not need the N-squared cell system.

::

    import espressomd.electrostatics
    mmm1d = espressomd.electrostatics.MMM1DGPU(prefactor=C, far_switch_radius=fr,
                                               maxPWerror=err, tune=False, bessel_cutoff=bc)
    mmm1d = espressomd.electrostatics.MMM1DGPU(prefactor=C, maxPWerror=err)

The first form sets parameters manually. The switch radius determines at which
xy-distance the force calculation switches from the near to the far
formula. If the Bessel cutoff is not explicitly given, it is determined
from the maximal pairwise error, otherwise this error only counts for
the near formula. The second tuning form just takes the maximal pairwise
error and tries out a lot of switching radii to find out the fastest one.

For details on the MMM family of algorithms, refer to appendix
:ref:`The MMM family of algorithms`.


.. _ScaFaCoS electrostatics:

ScaFaCoS electrostatics
-----------------------

:class:`espressomd.electrostatics.Scafacos`

|es| can use the methods from the ScaFaCoS *Scalable fast Coulomb solvers*
library. The specific methods available depend on the compile-time options of
the library, and can be queried using :meth:`espressomd.scafacos.available_methods`.

To use ScaFaCoS, create an instance of :class:`~espressomd.electrostatics.Scafacos`
and add it to the list of active actors. Three parameters have to be specified:
``prefactor`` (as defined in :eq:`coulomb_prefactor`), ``method_name``,
``method_params``. The method-specific parameters are described in the
ScaFaCoS manual. In addition, methods supporting tuning have a parameter
``tolerance_field`` which sets the desired root mean square accuracy for
the electric field.

To use a specific electrostatics solver from ScaFaCoS for your system,
e.g. ``ewald``, set its cutoff to :math:`1.5` and tune the other parameters
for an accuracy of :math:`10^{-3}`::

   import espressomd.electrostatics
   scafacos = espressomd.electrostatics.Scafacos(
      prefactor=1, method_name="ewald",
      method_params={"ewald_r_cut": 1.5, "tolerance_field": 1e-3})
   system.actors.add(scafacos)

For details of the various methods and their parameters please refer to
the ScaFaCoS manual. To use this feature, ScaFaCoS has to be built as a
shared library.
