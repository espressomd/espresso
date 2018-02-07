.. _Magnetostatics / Dipolar interactions:

Magnetostatics / Dipolar interactions
=====================================

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

where :math:`\mu_0` and :math:`\mu` are the vacuum permittivity and the relative permittivity of the background material, respectively.

Magnetostatic interactions are activated via the actor framework::

    from espressomd.magnetostatics import DipolarDirectSumCpu

    direct_sum=DipolarDirectSumCpu(prefactor=1)
    system.actors.add(direct_sum)
    #...
    system.actors.remove(direct_sum)

The magnetostatics algorithms for periodic boundary conditions require
some knowledge to use them properly. Uneducated use can result in
completely unphysical simulations.



.. _Dipolar P3M:

Dipolar P3M
~~~~~~~~~~~
This is the dipolar version of the P3M algorithm, described in :cite:`cerda08d`.
It is interfaced via :class:`espressomd.magnetostatics.DipolarP3M`.

Make sure that you know the relevance of the P3M parameters before using
P3M! If you are not sure, read the following references

Note that dipolar P3M does not work with non-cubic boxes.


The parameters of the dipolar P3M method can be tuned automatically, by providing `accuracy=<TARGET_ACCURACY>` to the method. 
It is also possible to pass a subset of the method parameters such as `mesh`. In that case, only the omitted parameters are tuned::


    import espressomd.magnetostatics as magnetostatics        
    p3m = magnetostatics.DipolarP3M(prefactor=1, mesh=32, accuracy=1E-4)
    system.actors.add(p3m)

It is important to note that the error estimates given in :cite:`cerda08a` used in the tuning contain assumptions about the system. In particular, a homogeneous system is assumed. If this is no longer the case during the simulation, actual force and torque errors can be significantly larger.

.. _Dipolar Layer Correction (DLC):


Dipolar Layer Correction (DLC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`espressomd.magnetostatic_extensions.DLC` 

The dipolar layer correction (DLC) is used in conjunction with the dipolar P3M method to calculate dipolar interactions in a 2D-periodic system.
It is based on :cite:`brodka04a` and the dipolar version of 
:ref:`Electrostatic Layer Correction (ELC)`.

Usage notes:

  * The non-periodic direction is always the `z`-direction.
  
  * The method relies on a slab of the simulation box perpendicular to the z-direction not to contain particles. The size in z-direction of this slab is controlled by the `gap_size` parameter. The user has to ensure that no particles enter this region by menas of constraints or by fixing the particles' z-coordinate. When there is no empty slab of the specified size, the method will silently produce wrong results.

  * The method can be tuned using the `accuracy` parameter. In contrast to the elctrostatic method, it refers to the energy. Furthermore, it is assumed that all dipole moment are as larger as the largest of the dipoles in the system. 

The method is used as follows::

    import espressomd.magnetostatics as magnetostatics
    import espressomd.magnetostatic_extensions as magnetostatic_extensions
    
    p3m = magnetostatics.DipolarP3M(prefactor=1, accuracy=1E-4)
    dlc = magnetostatic_extensions.DLC(maxPWerror=1E-5, gap_size=2.)
    system.actors.add(p3m)
    system.actors.add(dlc)



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


