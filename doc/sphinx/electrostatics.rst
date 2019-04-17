.. _Electrostatics:

Electrostatics
==============

The Coulomb (or electrostatic) interaction is defined as
follows. For a pair of particles at distance :math:`r` with charges
:math:`q_1` and :math:`q_2`, the interaction is given by

.. math:: U_C(r)=C \cdot \frac{q_1 q_2}{r}.

where

.. math::
   C=\frac{1}{4\pi \epsilon_0 \epsilon_r}
   :label: coulomb_prefactor

is a prefactor which can be set by the user. The commonly used Bjerrum length
:math:`l_B = e_o^2 / (4 \pi \epsilon_0 \epsilon_r k_B T)` is the length at
which the Coulomb energy between two unit charges is equal to the thermal
energy :math:`k_B T`.
Based on the this length, the prefactor is given by :math:`C=l_B k_B T`.

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

This example shows the general usage of an electrostatic method ``<SOLVER>``.
All of them need the Bjerrum length and a set of other required parameters.
First, an instance of the solver is created and only after adding it to the actors
list, it is activated. Internally the method calls a tuning routine on
activation to achieve the given accuracy::

    import espressomd
    from espressomd import electrostatics

    system = espressomd.System()
    solver = electrostatics.<SOLVER>(prefactor=C, <ADDITIONAL REQUIRED PARAMETERS>)
    system.actors.add(solver)

where the prefactor :math:`C` is defined as in Eqn. :eq:`coulomb_prefactor`

.. _Coulomb P3M:

Coulomb P3M
-----------

:class:`espressomd.electrostatics.P3M`

Required parameters:
    * ``prefactor``
    * ``accuracy``

For this feature to work, you need to have the ``fftw3`` library
installed on your system. In |es|, you can check if it is compiled in by
checking for the feature ``FFTW`` with ``espressomd.features``.
P3M requires full periodicity (1 1 1). Make sure that you know the relevance of the
P3M parameters before using P3M! If you are not sure, read the following
references
:cite:`ewald21,hockney88,kolafa92,deserno98,deserno98a,deserno00,deserno00a,cerda08a`.

.. _Tuning Coulomb P3M:

Tuning Coulomb P3M
~~~~~~~~~~~~~~~~~~

The tuning method is called when the handle of the Coulomb P3M is added to the
actor list. At this point, the system should already contain the charged
particles. Set parameters are fixed and not changed by the tuning algorithm.
This can be useful to speed up the tuning during testing or if the parameters
are already known.

To prevent the automatic tuning, set the ``tune`` parameter to ``False``.
To manually tune or retune P3M, call :meth:`espresso.electrostatics.P3M.Tune`.
Note, however, that this is a method the P3M object inherited from
:attr:`espressomd.electrostatics.ElectrostaticInteraction`.
All parameters passed to the method are fixed in the tuning routine. If not
specified in the ``Tune()`` method, the parameters ``prefactor`` and
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
the P3M method :cite:`hockney88` and its real space error :cite:`kolafa92` to
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

Required parameters:
    * ``prefactor``
    * ``accuracy``

The GPU implementation of P3M calculates the far field portion on the GPU.
It uses the same parameters and interface functionality as the CPU version of
the solver. It should be noted that this does not always provide significant
increase in performance. Furthermore it computes the far field interactions
with only single precision which limits the maximum precision. The algorithm
does not work in combination with the electrostatic extensions :ref:`Dielectric interfaces with the ICC algorithm`
and :ref:`Electrostatic Layer Correction (ELC)`.

.. _Debye-Hückel potential:

Debye-Hückel potential
----------------------

For a list of all parameters see :attr:`espressomd.electrostatics.DH`
Uses the Debye-Hückel electrostatic potential defined by

  .. math:: U^{C-DH} = C \cdot \frac{q_1 q_2 \exp(-\kappa r)}{r}\quad \mathrm{for}\quad r<r_{\mathrm{cut}}

where :math:`C` is defined as in Eqn. :eq:`coulomb_prefactor`.
The Debye-Hückel potential is an approximate method for calculating
electrostatic interactions, but technically it is treated as other
short-ranged non-bonding potentials. For :math:`r>r_{\mathrm cut}` it is
set to zero which introduces a step in energy. Therefore, it introduces
fluctuations in energy.

For :math:`\kappa = 0`, this corresponds to the plain Coulomb potential.


.. _Dielectric interfaces with the ICC algorithm:

Dielectric interfaces with the ICC\ :math:`\star` algorithm
-----------------------------------------------------------

The ICC\ :math:`\star` algorithm allows to take into account arbitrarily shaped
dielectric interfaces and dynamic charge induction. For instance, it can be
used to simulate a curved metallic boundary. This is done by iterating the
charge on a set of spatially fixed *ICC particles* until they correctly
represent the influence of the dielectric discontinuity. All *ICC particles*
need a certain area, normal vector and dielectric constant to specify the
surface. ICC relies on a Coulomb solver that is already initialized. So far, it
is implemented and well tested with the Coulomb solver P3M. ICC is an |es|
actor and can be activated via::

    icc = ICC(<See the following list of ICC parameters>)
    system.actors.add(icc)

Parameters are:

	* ``first_id``:
		ID of the first ICC Particle.
	* ``n_icc``:
		Total number of ICC Particles.
	* ``convergence``:
		Abort criteria of the iteration. It corresponds to the maximum relative
		change of any of the interface particle's charge.
	* ``relaxation``:
		SOR relaxation parameter.
	* ``ext_field``:
		Homogeneous electric field added to the calculation of dielectric boundary forces.
	* ``max_iterations``:
		Maximal number of iterations.
	* ``eps_out``:
		Relative permittivity of the outer region (where the particles are).
	* ``normals``:
		List of size ``n_icc`` with normal vectors pointing into the outer region.
	* ``areas``
		List of size ``n_icc`` with areas of the discretized surface.
	* ``sigmas``
		List of size ``n_icc`` with an additional surface charge density in
		absence of any charge induction
	* ``epsilons``
		List of size ``n_icc`` with the dielectric constant associated to the area.

The ICC particles are setup as normal |es| particles. Note that they should be
fixed in space and need an initial nonzero charge. The following usage example
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
            system.part.add(pos=[l * xi, l * yi, 0], q=-0.0001, fix=[1, 1, 1], type=icc_type)
    iccNormals.extend([0, 0, 1] * nicc_per_electrode)

    # Right electrode (normal [0,0,-1])
    for xi in xrange(nicc):
        for yi in xrange(nicc):
            system.part.add(pos=[l * xi, l * yi, box_l], q=0.0001, fix=[1, 1, 1], type=icc_type)
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
corresponding articles, mainly :cite:`espresso2,tyagi10a,kesselheim11a` before
using it.

.. _MMM2D:

MMM2D
-----

.. note::
    Required features: ``ELECTROSTATICS``, ``PARTIAL_PERIODIC``.

MMM2D is an electrostatics solver for explicit 2D periodic systems.
It can account for different dielectric jumps on both sides of the
non-periodic direction. MMM2D Coulomb method needs periodicity 1 1 0 and the
layered cell system. The performance of the method depends on the number of
slices of the cell system, which has to be tuned manually. It is
automatically ensured that the maximal pairwise error is smaller than
the given bound. Note that the user has to take care that the particles don't
leave the box in the non-periodic z-direction e.g. with constraints. By default,
no dielectric contrast is set and it is used as::

    mmm2d = electrostatics.MMM2D(prefactor=C, maxPWerror=1e-3)
    system.actors.add(mmm2d)

where the prefactor :math:`C` is defined in Eqn. :eq:`coulomb_prefactor`.
For a detailed list of parameters see :attr:`espressomd.electrostatics.MMM2D`.
The last two, mutually exclusive parameters ``dielectric`` and
``dielectric_constants_on`` allow to specify dielectric contrasts at the
upper and lower boundaries of the simulation box. The first form
specifies the respective dielectric constants in the media, which
however is only used to calculate the contrasts. That is, specifying
:math:`\epsilon_t=\epsilon_m=\epsilon_b=\text{const}` is always
identical to :math:`\epsilon_t=\epsilon_m=\epsilon_b=1`::

    mmm2d = electrostatics.MMM2D(prefactor=C, maxPWerror=1e-3, dielectric=1,
                                 top=1, mid=1, bot=1)

The second form specifies only the dielectric contrasts at the boundaries,
that is :math:`\Delta_t=\frac{\epsilon_m-\epsilon_t}{\epsilon_m+\epsilon_t}`
and :math:`\Delta_b=\frac{\epsilon_m-\epsilon_b}{\epsilon_m+\epsilon_b}`.
Using this form allows to choose :math:`\Delta_{t/b}=-1`, corresponding
to metallic boundary conditions::

    mmm2d = electrostatics.MMM2D(prefactor=C, maxPWerror=1e-3, dielectric_contrast_on=1,
                                 delta_mid_top=-1, delta_mid_bot=-1)

Using ``const_pot`` allows to maintain a constant electric potential difference ``pot_diff``
between the xy-planes at :math:`z=0` and :math:`z=L`, where :math:`L`
denotes the box length in :math:`z`-direction::

    mmm2d = electrostatics.MMM2D(prefactor=100.0, maxPWerror=1e-3, const_pot=1, pot_diff=100.0)

This is done by countering the total dipole moment of the system with the
electric field :math:`E_{induced}` and superposing a homogeneous electric field
:math:`E_{applied} = \frac{U}{L}` to retain :math:`U`. This mimics the
induction of surface charges :math:`\pm\sigma = E_{induced} \cdot \epsilon_0`
for planar electrodes at :math:`z=0` and :math:`z=L` in a capacitor connected
to a battery with voltage ``pot_diff``. Using 0 is equivalent to
:math:`\Delta_{t/b}=-1`.

Finally, the far cutoff setting should only be used for testing reasons,
otherwise you are more safe with the automatic tuning. If you even don't know
what it is, do not even think of touching the far cutoff. For details on the
MMM family of algorithms, refer to appendix :ref:`The MMM family of algorithms`.
Please cite :cite:`arnold02a` when using MMM2D.

A complete (but unphysical) sample script for a plate capacitor simulated with MMM2D
can be found in :file:`/samples/visualization_mmm2d.py`.

.. _Electrostatic Layer Correction (ELC):

Electrostatic Layer Correction (ELC)
------------------------------------

*ELC* can be used to simulate charged system with 2D periodicity. In more
detail, is a special procedure that converts a 3D electrostatic method to a 2D
method in computational order N. Currently, it only supports P3M. This means,
that you will first have to set up the P3M algorithm before using ELC. The
algorithm is definitely faster than MMM2D for larger numbers of particles
(:math:`>400` at reasonable accuracy requirements). The periodicity has to be
set to ``1 1 1`` still, *ELC* cancels the electrostatic contribution of the
periodic replica in **z-direction**. Make sure that you read the papers on ELC
(:cite:`arnold02c,icelc`) before using it. ELC is an |es| actor and is used
with::

    elc = electrostatic_extensions.ELC(gap_size=box_l * 0.2, maxPWerror=1e-3)
    system.actors.add(elc)


Parameters are:
    * ``gap_size``:
        The gap size gives the height of the empty region between the system box
        and the neighboring artificial images. |es| does not
        make sure that the gap is actually empty, this is the users
        responsibility. The method will compute fine if the condition is not
        fulfilled, however, the error bound will not be reached. Therefore you
        should really make sure that the gap region is empty (e.g. with wall
        constraints).
    * ``maxPWerror``:
        The maximal pairwise error sets the least upper bound (LUB) error of
        the force between any two charges without prefactors (see the papers).
        The algorithm tries to find parameters to meet this LUB requirements or
        will throw an error if there are none.
    * ``delta_mid_top``/``delta_mid_bot``:
        *ELC* can also be used to simulate 2D periodic systems with image charges,
        specified by dielectric contrasts on the non-periodic boundaries
        (:cite:`icelc`).  Similar to *MMM2D*, these can be set with the
        keywords ``delta_mid_bot`` and ``delta_mid_top``, setting the dielectric
        jump from the simulation region (*middle*) to *bottom* (at ``z<0``) and
        from *middle* to *top* (``z > box_l[2] - gap_size``). The fully metallic case
        ``delta_mid_top=delta_mid_bot=-1`` would lead to divergence of the
        forces/energies in *ELC* and is therefore only possible with the
        ``const_pot`` option.
    * ``const_pot``:
        As described, setting this to ``1`` leads to fully metallic boundaries and
        behaves just like the mmm2d parameter of the same name: It maintains a
        constant potential ``pot_diff`` by countering the total dipole moment of
        the system and adding a homogeneous electric field according to
        ``pot_diff``.
    * ``pot_diff``:
        Used in conjunction with ``const_pot`` set to 1, this sets the potential difference
        between the boundaries in the z-direction between ``z=0`` and
        ``z = box_l[2] - gap_size``.
    * ``far_cut``:
        The setting of the far cutoff is only intended for testing and allows to
        directly set the cutoff. In this case, the maximal pairwise error is
        ignored.
    * ``neutralize``:
        By default, ELC just as P3M adds a homogeneous neutralizing background
        to the system in case of a net charge. However, unlike in three dimensions,
        this background adds a parabolic potential across the
        slab :cite:`ballenegger09a`. Therefore, under normal circumstance, you will
        probably want to disable the neutralization for non-neutral systems.
        This corresponds then to a formal regularization of the forces and
        energies :cite:`ballenegger09a`. Also, if you add neutralizing walls
        explicitly as constraints, you have to disable the neutralization.
        When using a dielectric contrast or full metallic walls
        (``delta_mid_top != 0`` or ``delta_mid_bot != 0`` or
        ``const_pot=1``), ``neutralize`` is overwritten and switched off internally.
        Note that the special case of non-neutral systems with a *non-metallic* dielectric jump (eg.
        ``delta_mid_top`` or ``delta_mid_bot`` in ``]-1,1[``) is not covered by the
        algorithm and will throw an error.


.. _MMM1D:

MMM1D
-----

.. note::
    Required features: ``ELECTROSTATICS``, ``PARTIAL_PERIODIC`` for MMM1D, the GPU version additionally needs
    the features ``CUDA`` and ``MMM1D_GPU``.

::

    from espressomd.electrostatics import MMM1D
    from espressomd.electrostatics import MMM1DGPU

Please cite :cite:`arnold05a`  when using MMM1D.

See :attr:`espressomd.electrostatics.MMM1D` or
:attr:`espressomd.electrostatics.MMM1DGPU` for the list of available
parameters.

::

    mmm1d = MMM1D(prefactor=C, far_switch_radius=fr, maxPWerror=err, tune=False,
                  bessel_cutoff=bc)
    mmm1d = MMM1D(prefactor=C, maxPWerror=err)

where the prefactor :math:`C` is defined in Eqn. :eq:`coulomb_prefactor`.
MMM1D Coulomb method for systems with periodicity 0 0 1. Needs the
nsquared cell system (see section :ref:`Cellsystems`). The first form sets parameters
manually. The switch radius determines at which xy-distance the force
calculation switches from the near to the far formula. The Bessel cutoff
does not need to be specified as it is automatically determined from the
particle distances and maximal pairwise error. The second tuning form
just takes the maximal pairwise error and tries out a lot of switching
radii to find out the fastest one. If this takes too long, you can
change the value of the setmd variable ``timings``, which controls the number of
test force calculations.

::

    mmm1d_gpu = MMM1DGPU(prefactor=C, far_switch_radius=fr, maxPWerror=err,
                         tune=False, bessel_cutoff=bc)
    mmm1d_gpu = MMM1DGPU(prefactor=C, maxPWerror=err)

MMM1D is also available in a GPU implementation. Unlike its CPU
counterpart, it does not need the nsquared cell system. The first form
sets parameters manually. The switch radius determines at which
xy-distance the force calculation switches from the near to the far
formula. If the Bessel cutoff is not explicitly given, it is determined
from the maximal pairwise error, otherwise this error only counts for
the near formula. The second tuning form just takes the maximal pairwise
error and tries out a lot of switching radii to find out the fastest
one.

For details on the MMM family of algorithms, refer to appendix :ref:`The MMM family of algorithms`.


.. _Scafacos Electrostatics:

Scafacos Electrostatics
-----------------------

Espresso can use the electrostatics methods from the SCAFACOS *Scalable
fast Coulomb solvers* library. The specific methods available depend on the compile-time options of the library, and can be queried using :attr:`espressomd.scafacos.available_methods`

To use SCAFACOS, create an instance of :attr:`espressomd.electrostatics.Scafacos` and add it to the list of active actors. Three parameters have to be specified:

* ``method_name``: name of the SCAFACOS method being used.
* ``method_params``: dictionary containing the method-specific parameters
* ``prefactor``: Coulomb prefactor as defined in :eq:`coulomb_prefactor`.

The method-specific parameters are described in the SCAFACOS manual.
Additionally, methods supporting tuning have the parameter ``tolerance_field`` which sets the desired root mean square accuracy for the electric field

To use the, e.g.,  ``ewald`` solver from SCAFACOS as electrostatics solver for your system, set its
cutoff to :math:`1.5` and tune the other parameters for an accuracy of
:math:`10^{-3}`, use::

  from espressomd.electrostatics import Scafacos
  scafacos = Scafacos(prefactor=1, method_name="ewald",
                      method_params={"ewald_r_cut": 1.5, "tolerance_field": 1e-3})
  system.actors.add(scafacos)


For details of the various methods and their parameters please refer to
the SCAFACOS manual. To use this feature, SCAFACOS has to be built as a shared library. SCAFACOS can be used only once, either for Coulomb or for dipolar interactions.

