.. _Magnetostatics / Dipolar interactions:

Magnetostatics / Dipolar interactions
=====================================

.. _Dipolar interaction:

Dipolar interaction
-------------------

inter magnetic 0.0 inter magnetic inter magnetic

These commands can be used to set up magnetostatic interactions, which
is defined as follows:

.. math::

   U^{D-P3M}(\vec{r}) = D \cdot \left( \frac{(\vec{\mu}_i \cdot \vec{\mu}_j)}{r^3} 
     - \frac{3  (\vec{\mu}_i \cdot \vec{r})  (\vec{\mu}_j \cdot \vec{r}) }{r^5} \right)

where :math:`r=|\vec{r}|`.
The prefactor :math:`D` is can be set by the user and is given by

.. math::

  D =\frac{\mu_0 \mu}{4\pi}

Computing magnetostatic interactions is computationally very expensive.
features some state-of-the-art algorithms to deal with these
interactions as efficiently as possible, but almost all of them require
some knowledge to use them properly. Uneducated use can result in
completely unphysical simulations.

The commands above work as their counterparts for the electrostatic
interactions (see section ). Variant disables dipolar interactions.
Variant returns the current parameters of the dipolar interaction as a
Tcl-list using the same syntax as used to setup the method,

coulomb 1.0 p3m 7.75 8 5 0.1138 0.0 coulomb epsilon 0.1 n_interpol
32768 mesh_off 0.5 0.5 0.5

Variant is the generic syntax to set up a specific method or its
parameters, the details of which are described in the following
subsections. Note that using the magnetostatic interaction also requires
assigning dipole moments to the particles. This is done using the
``part`` command to set the dipole moment ``dip``,

inter coulomb 1.0 p3m tune accuracy 1e-4 part 0 dip 1 0 0; part 1 dip 0
0 1

.. _Dipolar P3M:

Dipolar P3M
~~~~~~~~~~~

inter magnetic p3m

This command activates the P3M method to compute the dipolar
interactions between charged particles. The different parameters are
described in more detail in :cite:`cerda08a`.

    The real space cutoff as a positive floating point number.

    The number of mesh points, as a single positive integer.

    The *charge-assignment order*, an integer between :math:`0` and
    :math:`7`.

    The Ewald parameter as a positive floating point number.

Make sure that you know the relevance of the P3M parameters before using
P3M! If you are not sure, read the following references
:cite:`ewald21,hockney88,kolafa92,deserno98,deserno98a,deserno00,deserno00a`.

Note that dipolar P3M does not work with non-cubic boxes.

.. _Tuning dipolar P3M:

Tuning dipolar P3M
^^^^^^^^^^^^^^^^^^

| inter magnetic p3m accuracy

Tuning dipolar P3M works exactly as tuning Coulomb P3M. Therefore, for
details on how to tune the algorithm, refer to the documentation of
Coulomb P3M (see section ).

For the magnetic case, the expressions of the error estimate are given
in :cite:`cerda08a`.

.. _Dipolar Layer Correction (DLC):

Dipolar Layer Correction (DLC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inter magnetic mdlc

Like ELC but applied to the case of magnetic dipoles, but here the
accuracy is the one you wish for computing the energy. is set to a value
that, assuming all dipoles to be as larger as the largest of the dipoles
in the system, the error for the energy would be smaller than the value
given by accuracy. At this moment you cannot compute the accuracy for
the forces, or torques, nonetheless, usually you will have an error for
forces and torques smaller than for energies. Thus, the error for the
energies is an upper boundary to all errors in the calculations.

At present, the program assumes that the gap without particles is along
the z-direction. The gap-size is the length along the z-direction of the
volume where particles are not allowed to enter.

As a reference for the DLC method, see :cite:`brodka04a`.

.. _Dipolar all-with-all and no replicas (DAWAANR):

Dipolar all-with-all and no replicas (DAWAANR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inter magnetic dawaanr

This interaction calculates energies and forces between dipoles by
explicitly summing over all pairs. For the directions in which the
system is periodic (as defined by ``setmd periodic``), it applies the
minimum image convention, i.e. the interaction is effectively cut off at
half a box length.

In periodic systems, this method should only be used if it is not
possible to use dipolar P3M or DLC, because those methods have a far
better accuracy and are much faster. In a non-periodic system, the
DAWAANR-method gives the exact result.

.. _Magnetic Dipolar Direct Sum (MDDS) on CPU:

Magnetic Dipolar Direct Sum (MDDS) on CPU
-----------------------------------------

inter magnetic mdds n_cut

The command enables the “magnetic dipolar direct sum”. The dipole-dipole
interaction is computed by explicitly summing over all pairs. If the
system is periodic in one or more directions, the interactions with
further replicas of the system in all periodic directions is explicitly
computed.

As it is very slow, this method is not intended to do simulations, but
rather to check the results you get from more efficient methods like
P3M.

.. _Dipolar direct sum on gpu:

Dipolar direct sum on gpu
-------------------------

This interaction calculates energies and forces between dipoles by
explicitly summing over all pairs. For the directions in which the
system is periodic (as defined by ``setmd periodic``), it applies the
minimum image convention, i.e. the interaction is effectively cut off at
half a box length.

The calculations are performed on the gpu in single precision. The
implementation is optimized for large systems of several thousand
particles. It makes use of one thread per particle. When there are fewer
particles than the number of threads the gpu can execute simultaneously,
the rest of the gpu remains idle. Hence, the method will perform poorly
for small systems.

To use the method, create an instance of :attr:`espressomd.magnetostatics.DipolarDirectSumGpu` and add it to the system's list of active actors. The only required parameter is the Bjerrum length::
  
  from espressomd.magnetostatics import DipolarDirectSumGpu
  dds=DipolarDirectSumGpu(bjerrum_length=1)
  system.actors.add(dds)

.. _Barnes-Hut octree sum on gpu:

Barnes-Hut octree sum on gpu
----------------------------

This interaction calculates energies and forces between dipoles by
summing over the spatial octree cells (aka ``leaves``).
Far enough cells are considered as a single dipole with a cumulative
vector in the cell center of mass. Parameters which determine that the
cell is far enough are :math:`I_{\mathrm{tol}}^2` and
:math:`\varepsilon^2` which define a fraction of the cell and
an additive distance respectively. For the detailed description of the
Barnes-Hut method application to the dipole-dipole interactions, please
refer to :cite:`Polyakov2013`.

To use the method, create an instance of :attr:`espressomd.magnetostatics.DipolarBarnesHutGpu` and add it to the system's list of active actors::
  
  from espressomd.magnetostatics import DipolarBarnesHutGpu
  bh=DipolarBarnesHutGpu(prefactor = pf_dds_gpu, epssq = 200.0, itolsq = 8.0)
  system.actors.add(bh)

.. _Scafacos Magnetostatics:

Scafacos Magnetostatics
-----------------------

Espresso can use the methods from the Scafacos *Scalable fast Coulomb
solvers* library for dipoles, if the methods support dipolar
calculations. The feature SCAFACOS_DIPOLES has to be added to
myconfig.hpp to activate this feature. At the time of this writing (May
2017) dipolar calculations are only included in the ``dipolar`` branch of the Scafacos code.

To use SCAFACOS, create an instance of :attr:`espressomd.magnetostatics.Scafacos` and add it to the list of active actors. Three parameters have to be specified:
* method_name: name of the SCAFACOS method being used.
* method_params: dictionary containing the method-specific parameters
* bjerrum_length
The method-specific parameters are described in the SCAFACOS manual.
Additionally, methods supporting tuning have the parameter ``tolerance_field`` which sets the desired root mean square accuracy for the electric field 

For details of the various methods and their parameters please refer to
the SCAFACOS manual. To use this feature, SCAFACOS has to be built as a shared library. SCAFACOS can be used only once, either for coulomb or for dipolar interactions.


