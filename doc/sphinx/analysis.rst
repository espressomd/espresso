.. _Analysis:

Analysis
========

|es| provides two concepts of system analysis:

- :ref:`Direct analysis routines`: The :mod:`espressomd.analyze` module provides
  online-calculation of specialized local and global observables with
  calculation and data accumulation performed in the core.
- :ref:`Observables and correlators`: This provides a more flexible concept of
  in-core analysis, where a certain observable (:ref:`Available observables`),
  a rule for data accumulation (:ref:`Accumulators`) and/or correlation (:ref:`Correlations`) can be defined.


.. _Direct analysis routines:

Direct analysis routines
------------------------

The direct analysis commands can be classified into two types:

- Instantaneous analysis routines, that only take into account the current configuration of the system:

    - :ref:`Energies`
    - :ref:`Pressure`
    - :ref:`Momentum of the System`
    - :ref:`Minimal distances between particles`
    - :ref:`Particles in the neighborhood`
    - :ref:`Particle distribution`
    - :ref:`Radial distribution function` with ``rdf_type='rdf'``
    - :ref:`Structure factor`
    - :ref:`Center of mass`
    - :ref:`Moment of inertia matrix`
    - :ref:`Gyration tensor`
    - :ref:`Stress Tensor`

- Analysis on stored configurations, added by :meth:`espressomd.analyze.Analysis.append`:
    - :ref:`Radial distribution function` with ``rdf_type='<rdf>'``
    - :ref:`Chains`

.. _Energies:

Energies
~~~~~~~~
:meth:`espressomd.analyze.Analysis.energy`

Returns the energies of the system.
The different energetic contributions to the total energy can also be obtained (kinetic, bonded,non-bonded, Coulomb).

For example, ::

>>> energy = system.analysis.energy()
>>> print(energy["total"])
>>> print(energy["kinetic"])
>>> print(energy["bonded"])
>>> print(energy["non_bonded"])


.. _Momentum of the system:

Momentum of the System
~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.linear_momentum`

This command returns the total linear momentum of the particles and the
lattice-Boltzmann (LB) fluid, if one exists. Giving the optional
parameters either causes the command to ignore the contribution of LB or
of the particles.

.. _Minimal distances between particles:

Minimal distances between particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`espressomd.analyze.Analysis.min_dist`
Returns the minimal distance between all particles in the system.

When used with type-lists as arguments, then the minimal distance between particles of only those types is determined.


For example, ::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[100, 100, 100])
    >>> for i in range(10):
    ...     system.part.add(id=i, pos=[1.0, 1.0, i**2], type=0)
    >>> system.analysis.min_dist()
    1.0


.. _Particles in the neighborhood:

Particles in the neighborhood
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`espressomd.analyze.Analysis.nbhood`

Returns a list of the particle ids of that fall within a given radius of a target position.
For example, ::

    idlist = system.analysis.nbhood(pos=system.box_l * 0.5, r_catch=5.0)

.. _Particle distribution:

Particle distribution
~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.distribution`

Returns the distance distribution of particles
(probability of finding a particle of a certain type at a specified distance around
a particle of another specified type, disregarding the fact that a spherical shell of a
larger radius covers a larger volume).
The distance is defined as the *minimal* distance between a particle of one group to any of the other
group.

Two arrays are returned corresponding to the normalized distribution and the bins midpoints, for example ::

    >>> system = espressomd.System(box_l=[10, 10, 10])
    >>> for i in range(5):
    ...     system.part.add(id=i, pos=i * system.box_l, type=0)
    >>> bins, count = system.analysis.distribution(type_list_a=[0], type_list_b=[0],
    ...                                            r_min=0.0, r_max=10.0, r_bins=10)
    >>>
    >>> print(bins)
    [ 0.5  1.5  2.5  3.5  4.5  5.5  6.5  7.5  8.5  9.5]
    >>> print(count)
    [ 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]


.. _Radial distribution function:

Radial distribution function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.rdf`

Calculates a radial distribution function for given particle type and binning.
The ``rdf_type`` defines if the analysis is performed on the current configuration (``rdf_type='rdf'``)
or on averaged configurations stored with :meth:`analyze.append() <espressomd.analyze.Analysis.append>` (``rdf_type='<rdf>'``).

For example, ::

    rdf_bins = 100
    r_min = 0.0
    r_max = system.box_l[0] / 2.0
    r, rdf_01 = S.analysis.rdf(rdf_type='<rdf>', type_list_a=[0], type_list_b=[1],
                               r_min=r_min, r_max=r_max, r_bins=rdf_bins)
    rdf_fp = open("rdf.dat", 'w')
    for i in range(rdf_bins):
        rdf_fp.write("%1.5e %1.5e %1.5e %1.5e\n" % (r[i], rdf_01[i]))
    rdf_fp.close()


.. _Structure factor:

Structure factor
~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.structure_factor`

Calculate the structure factor for given types.

Returns the spherically averaged structure factor :math:`S(q)` of
particles specified in ``sf_types``. :math:`S(q)` is calculated for all possible
wave vectors :math:`\frac{2\pi}{L} \leq q \leq \frac{2\pi}{L}` up to ``sf_order``.


.. _Center of mass:

Center of mass
~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.center_of_mass`

Returns the center of mass of particles of the given type given by ``part_type``.


.. _Moment of inertia matrix:

Moment of inertia matrix
~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.moment_of_inertia_matrix`

Returns the 3x3 moment of inertia matrix for particles of a given type.


.. _Gyration tensor:

Gyration tensor
~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.gyration_tensor`

Analyze the gyration tensor of particles of a given type, or of all particles in the system if no type is given. Returns a dictionary containing the squared radius of gyration, three shape descriptors (asphericity, acylindricity, and relative shape anisotropy), eigenvalues of the gyration tensor and their corresponding eigenvectors. The eigenvalues are sorted in descending order.


.. _Pressure:

Pressure
~~~~~~~~

:meth:`espressomd.analyze.Analysis.pressure`

Computes the instantaneous virial pressure for an isotropic and homogeneous system. It
returns all the contributions to the total pressure as well as the total pressure (see :meth:`espressomd.analyze.Analysis.pressure`).

The instantaneous pressure is calculated (if there are no electrostatic interactions)
by the volume averaged, direction averaged instantaneous virial pressure

.. math::
     p = \frac{2E_{\text{kinetic}}}{Vf} + \frac{\sum_{j>i} {F_{ij}r_{ij}}}{3V}
     :label: eqptens

where :math:`f=3` is the number of translational degrees of freedom of
each particle, :math:`V` is the volume of the system,
:math:`E_{\text{kinetic}}` is the kinetic energy, :math:`F_{ij}` the force
between particles i and j, and :math:`r_{ij}` is the distance between
them. The kinetic energy divided by the degrees of freedom is

.. math:: \frac{2E_{\text{kinetic}}}{f} = \frac{1}{3}\sum_{i} {m_{i}v_{i}^{2}}.

Note that Equation :eq:`eqptens` can only be applied to pair potentials and
central forces. Description of how contributions from other interactions
are calculated is beyond the scope of this manual. Three body potentials
are implemented following the procedure in
Ref. :cite:`thompson09a`. A different formula is used to
calculate contribution from electrostatic interactions. For
electrostatic interactions in P3M, the :math:`k`-space contribution is implemented according to :cite:`essmann95a`.
The implementation of the Coulomb P3M pressure is tested against LAMMPS.

Four-body dihedral potentials are not included. Except of
``VIRTUAL_SITES_RELATIVE`` constraints all other
constraints of any kind are not currently accounted for in the pressure
calculations. The pressure is no longer correct, e.g., when particles
are confined to a plane.

Note: The different contributions which are returned are the summands that arise from force splitting :math:`\vec{F}_{i,j}={\vec{F}_{i,j}}_\text{bonded}+{\vec{F}_{i,j}}_\text{nonbonded}+...` in the virial pressure formula. Later when the user calculates the ensemble average via e.g. :math:`\langle p \rangle \approx 1/N \sum_{i=1}^N p_i` however the ensemble average with all interactions present is performed. That means the contributions are not easy to interpret! Those are the contributions to the stress/pressure in a system where all interactions are present and therefore in a coupled system.

.. _Stress Tensor:

Stress Tensor
~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.stress_tensor`

Computes the volume averaged instantaneous stress tensor of the system with options which are
described by in :meth:`espressomd.analyze.Analysis.stress_tensor`. It is called a stress tensor but the sign convention follows that of a pressure tensor.
In general do only use it for (on average) homogeneous systems. For inhomogeneous systems you need to use the local stress tensor.

The instantaneous virial stress tensor is calculated by

.. math:: p_{(k,l)} = \frac{\sum_{i} {m_{i}v_{i}^{(k)}v_{i}^{(l)}}}{V} + \frac{\sum_{j>i}{F_{ij}^{(k)}r_{ij}^{(l)}}}{V}

where the notation is the same as for the pressure. The superscripts :math:`k`
and :math:`l` correspond to the components in the tensors and vectors.

If electrostatic interactions are present then also the coulombic parts of the stress tensor need to be calculated. If P3M is present, then the instantaneous stress tensor is added to the above equation in accordance with :cite:`essmann95a` :

.. math :: p^\text{Coulomb, P3M}_{(k,l)} =p^\text{Coulomb, P3M, dir}_{(k,l)} + p^\text{Coulomb, P3M, rec}_{(k,l)},

where the first summand is the short ranged part and the second summand is the long ranged part.

The short ranged part is given by:

.. math :: p^\text{Coulomb, P3M, dir}_{(k,l)}= \frac{1}{4\pi \varepsilon_0 \varepsilon_r} \frac{1}{2V} \sum_{\vec{n}}^* \sum_{i,j=1}^N q_i q_j \left( \frac{ \mathrm{erfc}(\beta |\vec{r}_j-\vec{r}_i+\vec{n}|)}{|\vec{r}_j-\vec{r}_i+\vec{n}|^3} + \\ \frac{2\beta \pi^{-1/2} \exp(-(\beta |\vec{r}_j-\vec{r}_i+\vec{n}|)^2)}{|\vec{r}_j-\vec{r}_i+\vec{n}|^2} \right) (\vec{r}_j-\vec{r}_i+\vec{n})_k (\vec{r}_j-\vec{r}_i+\vec{n})_l,

where :math:`\beta` is the P3M splitting parameter, :math:`\vec{n}` identifies the periodic images, the asterisk denotes that terms with :math:`\vec{n}=\vec{0}` and i=j are omitted.
The long ranged (k-space) part is given by:

.. math :: p^\text{Coulomb, P3M, rec}_{(k,l)}= \frac{1}{4\pi \varepsilon_0 \varepsilon_r} \frac{1}{2 \pi V^2} \sum_{\vec{k} \neq \vec{0}} \frac{\exp(-\pi^2 \vec{k}^2/\beta^2)}{\vec{k}^2} |S(\vec{k})|^2 \cdot (\delta_{k,l}-2\frac{1+\pi^2\vec{k}^2/\beta^2}{\vec{k}^2} \vec{k}_k \vec{k}_l),

where :math:`S(\vec{k})` is the Fourier transformed charge density. Compared to Essmann we do not have the contribution :math:`p^\text{corr}_{k,l}` since we want to calculate the pressure that arises from all particles in the system.

Note: The different contributions which are returned are the summands that arise from force splitting :math:`\vec{F}_{i,j}={\vec{F}_{i,j}}_\text{bonded}+{\vec{F}_{i,j}}_\text{nonbonded}+...` in the virial stress tensor formula.
Later when the user calculates the stress tensor via :math:`\langle p_{(k,l)}\rangle  \approx 1/N \sum_{i=1}^N p_{k,l}` however the ensemble average with all interactions present is performed.
That means the contributions are not easy to interpret! Those are the contributions to the stress/pressure in a system where all interactions are present and therefore in a coupled system.

Note that the angular velocities of the particles are not included in
the calculation of the stress tensor.

.. _Chains:

Chains
~~~~~~

All analysis functions in this section require the topology of the chains to be set correctly.
The above set of functions is designed to facilitate analysis of molecules.
Molecules are expected to be a group of particles comprising a contiguous range of particle IDs.
Each molecule is a set of consecutively numbered particles and all molecules are supposed to consist of the same number of particles.

Some functions in this group require that the particles constituting a molecule are connected into
linear chains (particle :math:`n` is connected to :math:`n+1` and so on)
while others are applicable to molecules of whatever topology.


.. _End to end distance:

End-to-end distance
^^^^^^^^^^^^^^^^^^^
:meth:`espressomd.analyze.Analysis.calc_re`

Returns the quadratic end-to-end-distance and its root averaged over all chains.

.. _Radius of gyration:

Radius of gyration
^^^^^^^^^^^^^^^^^^
:meth:`espressomd.analyze.Analysis.calc_rg`

Returns the radius of gyration averaged over all chains.
It is a radius of a sphere, which would have the same moment of inertia as the
molecule, defined as

.. math::

   \label{eq:Rg}
   R_{\mathrm G}^2 = \frac{1}{N} \sum\limits_{i=1}^{N} \left(\vec r_i - \vec r_{\mathrm{cm}}\right)^2\,,

where :math:`\vec r_i` are position vectors of individual particles
constituting a molecule and :math:`\vec r_{\mathrm{cm}}` is the position
vector of its center of mass. The sum runs over all :math:`N` particles
comprising the molecule. For more information see any polymer science
book, e.g. :cite:`rubinstein03a`.


.. _Hydrodynamic radius:

Hydrodynamic radius
^^^^^^^^^^^^^^^^^^^
:meth:`espressomd.analyze.Analysis.calc_rh`

Returns the hydrodynamic radius averaged over all chains.
The following formula is used for the computation:

.. math::

   \label{eq:Rh}
   \frac{1}{R_{\mathrm H}} = \frac{2}{N(N-1)} \sum\limits_{i=1}^{N} \sum\limits_{j<i}^{N} \frac{1}{|\vec r_i - \vec r_j|}\,,

The above-mentioned formula is only valid under certain assumptions. For
more information, see Chapter 4 and equation 4.102
in :cite:`doi86a`.
Note that the hydrodynamic radius is sometimes defined in a similar fashion but with a denominator of :math:`N^2` instead of :math:`N(N-1)` in the prefactor.
Both versions are equivalent in the :math:`N\rightarrow \infty` limit but give numerically different values for finite polymers.


.. _Observables and correlators:

Observables and correlators
---------------------------

Observables extract properties of the particles and the LB fluid and
return either the raw data or a statistic derived from them. The
Observable framework is progressively replacing the Analysis framework.
This is motivated by the fact, that sometimes it is desirable that the
analysis functions do more than just return a value to the scripting
interface. For some observables it is desirable to be sampled every few
integration steps. In addition, it should be possible to pass the
observable values to other functions which compute history-dependent
quantities, such as correlation functions. All this should be done
without the need to interrupt the integration by passing the control to
the script level and back, which produces a significant overhead when
performed too often.

Some observables in the core have their corresponding counterparts in
the :mod:`espressomd.analyze` module. However, only the core-observables can be used
on the fly with the toolbox of the correlator and on the fly analysis of
time series.
Similarly, some special cases of using the correlator have
their redundant counterparts in :mod:`espressomd.analyze`,
but the correlator provides a general and
versatile toolbox which can be used with any implemented
core-observables.

The first step of the core analysis is to create an observable.
An observable in the sense of the core analysis can be considered as a
rule how to compute a certain set of numbers from a given state of the
system or a role how to collect data from other observables. Any
observable is represented as a single array of double values in the core.
Any more complex shape (tensor, complex number, …) must be compatible to this
prerequisite. Every observable however documents the storage order and returns
a reshaped numpy array.

The observables can be used in parallel simulations. However,
not all observables carry out their calculations in parallel.
Instead, the entire particle configuration is collected on the head node, and the calculations are carried out there.
This is only performance-relevant if the number of processor cores is large and/or interactions are calculated very frequently.

.. _Using observables:

Using observables
~~~~~~~~~~~~~~~~~

The observables are represented as Python classes derived from
:class:`espressomd.observables.Observable`. They are contained in
the ``espressomd.observables`` module. An observable is instantiated as
follows

::

    from espressomd.observables import ParticlePositions
    part_pos = ParticlePositions(ids=(1, 2, 3, 4, 5))

Here, the keyword argument ``ids`` specifies the ids of the particles,
which the observable should take into account.

The current value of an observable can be obtained using its
:meth:`~espressomd.observables.Observable.calculate` method::

    print(part_pos.calculate())

Profile observables have additional methods
:meth:`~espressomd.observables.ProfileObservable.bin_centers` and
:meth:`~espressomd.observables.ProfileObservable.bin_edges` to facilitate
plotting of histogram slices with functions that require either bin centers
or bin edges for the axes. Example::

    import matplotlib.pyplot as plt
    import numpy as np
    import espressomd
    import espressomd.observables

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.part.add(id=0, pos=[4.0, 3.0, 6.0])
    system.part.add(id=1, pos=[7.0, 3.0, 6.0])

    # histogram in Cartesian coordinates
    density_profile = espressomd.observables.DensityProfile(
        ids=[0, 1],
        n_x_bins=8, min_x=1.0, max_x=9.0,
        n_y_bins=8, min_y=1.0, max_y=9.0,
        n_z_bins=4, min_z=4.0, max_z=8.0)
    obs_data = density_profile.calculate()
    obs_bins = density_profile.bin_centers()

    # 1D slice: requires bin centers
    plt.plot(obs_bins[:, 2, 2, 0], obs_data[:, 2, 2])
    plt.show()

    # 2D slice: requires extent
    plt.imshow(obs_data[:, :, 2].T, origin='lower',
               extent=[density_profile.min_x, density_profile.max_x,
                       density_profile.min_y, density_profile.max_y])
    plt.show()

    # histogram in cylindrical coordinates
    density_profile = espressomd.observables.CylindricalDensityProfile(
        ids=[0, 1], center=[5.0, 5.0, 0.0], axis=[0, 0, 1],
        n_r_bins=8, min_r=0.0, max_r=4.0,
        n_phi_bins=16, min_phi=-np.pi, max_phi=np.pi,
        n_z_bins=4, min_z=4.0, max_z=8.0)
    obs_data = density_profile.calculate()
    obs_bins = density_profile.bin_edges()

    # 2D slice: requires bin edges
    fig = plt.figure()
    ax = fig.add_subplot(111, polar='True')
    r = obs_bins[:, 0, 0, 0]
    phi = obs_bins[0, :, 0, 1]
    ax.pcolormesh(phi, r, obs_data[:, :, 2])
    plt.show()


.. _Available observables:

Available observables
~~~~~~~~~~~~~~~~~~~~~

The following list contains some of the available observables. You can find
documentation for all available observables in :mod:`espressomd.observables`.

- Observables working on a given set of particles:

   - :class:`~espressomd.observables.ParticlePositions`: Positions of the particles

   - :class:`~espressomd.observables.ParticleVelocities`: Velocities of the particles

   - :class:`~espressomd.observables.ParticleForces`: Forces on the particles

   - :class:`~espressomd.observables.ParticleBodyVelocities`: The particles' velocities in their respective body-fixed frames (as per their orientation in space stored in their quaternions).

   - :class:`~espressomd.observables.ParticleAngularVelocities`: The particles' angular velocities in the space-fixed frame

   - :class:`~espressomd.observables.ParticleBodyAngularVelocities`: As above, but in the particles' body-fixed frame.

- Observables working on a given set of particles and returning reduced quantities:

   - :class:`~espressomd.observables.Current`: Total current of the system

   - :class:`~espressomd.observables.DipoleMoment`: Total electric dipole moment of the system obtained based on unfolded positions

   - :class:`~espressomd.observables.MagneticDipoleMoment`: Total magnetic dipole moment of the system based on the :attr:`espressomd.particle_data.ParticleHandle.dip` property.

   - :class:`~espressomd.observables.ComPosition`: The system's center of mass based on unfolded coordinates

   - :class:`~espressomd.observables.ComVelocity`: Velocity of the center of mass

   - :class:`~espressomd.observables.ParticleDistances`: Distances between particles on a polymer chain.

   - :class:`~espressomd.observables.TotalForce`: Sum of the forces on the particles

   - :class:`~espressomd.observables.BondAngles`: Angles between bonds on a polymer chain.

   - :class:`~espressomd.observables.BondDihedrals`: Dihedral angles between bond triples on a polymer chain.

   - :class:`~espressomd.observables.CosPersistenceAngles`: Cosine of angles between bonds. The ``i``-th value in the result vector corresponds to the cosine of the angle between
     bonds that are separated by ``i`` bonds. This observable might be useful for measuring the persistence length of a polymer.

- Profile observables sampling the spatial profile of various quantities:

   - :class:`~espressomd.observables.DensityProfile`

   - :class:`~espressomd.observables.FluxDensityProfile`

   - :class:`~espressomd.observables.ForceDensityProfile`

   - :class:`~espressomd.observables.LBVelocityProfile`

- Observables sampling the cylindrical profile of various quantities:

   - :class:`~espressomd.observables.CylindricalDensityProfile`

   - :class:`~espressomd.observables.CylindricalFluxDensityProfile`

   - :class:`~espressomd.observables.CylindricalVelocityProfile`

   - :class:`~espressomd.observables.CylindricalLBFluxDensityProfileAtParticlePositions`

   - :class:`~espressomd.observables.CylindricalLBVelocityProfileAtParticlePositions`

- System-wide observables

   - :class:`~espressomd.observables.StressTensor`: Total stress tensor (see :ref:`stress tensor`)

   - :class:`~espressomd.observables.DPDStress`


.. _Correlations:

Correlations
~~~~~~~~~~~~

Time correlation functions are ubiquitous in statistical mechanics and
molecular simulations when dynamical properties of many-body systems are
concerned. A prominent example is the velocity autocorrelation function,
:math:`\left< \mathbf{v}(t) \cdot \mathbf{v}(t+\tau) \right>` which is
used in the Green-Kubo relations. In general, time correlation functions
are of the form

.. math::

   C(\tau) = \left<A\left(t\right) \otimes B\left(t+\tau\right)\right>


where :math:`t` is time, :math:`\tau` is the lag time (time difference)
between the measurements of (vector) observables :math:`A` and
:math:`B`, and :math:`\otimes` is an operator which produces the vector
quantity :math:`C` from :math:`A` and :math:`B`. The ensemble average
:math:`\left< \cdot \right>` is taken over all time origins \ :math:`t`.
Correlation functions describing dynamics of large and complex molecules
such as polymers span many orders of magnitude, ranging from MD time
step up to the total simulation time.

A correlator takes one or two observables, obtains values from them during the simulation and
finally uses a fast correlation algorithm which enables efficient computation
of correlation functions spanning many orders of magnitude in the lag time.

The implementation for computing averages and error estimates of a time series
of observables relies on estimates of autocorrelation functions and the
respective autocorrelation times. The correlator provides the same
functionality as a by-product of computing the correlation function.

An example of the usage of observables and correlations is provided in
the script :file:`samples/observables_correlators.py`.

.. _Creating a correlation:

Creating a correlation
^^^^^^^^^^^^^^^^^^^^^^

Each correlator is represented by an instance of the :class:`espressomd.accumulators.Correlator`. Please see its documentation for an explanation of the arguments that have to be passed to the constructor.

Correlators can be registered for automatic updating during the
integration by adding them to :attr:`espressomd.system.System.auto_update_accumulators`.

::

    system.auto_update_accumulators.add(corr)

Alternatively, an update can triggered by calling the ``update()`` method of the correlator instance. In that case, one has to make sure to call the update in the correct time intervals.


The current on-the-fly correlation result can of a correlator can be obtained using its ``result()`` method.
The final result (including the latest data in the buffers) is obtained using the ``finalize()`` method. After this, no further update of the correlator is possible.

.. _Example\: Calculating a particle's diffusion coefficient:

Example: Calculating a particle's diffusion coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For setting up an observable and correlator to obtain the mean square displacement of particle 0, use::

    pos_obs = ParticlePositions(ids=(0,))
    c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., delta_N=10,
                       corr_operation="square_distance_componentwise", compress1="discard1")

To obtain the velocity auto-correlation function of particle 0, use::

    obs = ParticleVelocities(ids=(0,))
    c_vel = Correlator(obs1=vel_obs, tau_lin=16, tau_max=20., delta_N=1,
                       corr_operation="scalar_product", compress1="discard1")

The full example can be found in :file:`samples/diffusion_coefficient.py`.


.. _Details of the multiple tau correlation algorithm:

Details of the multiple tau correlation algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we briefly describe the multiple tau correlator which is
implemented in |es|. For a more detailed description and discussion of its
behavior with respect to statistical and systematic errors, please read
the cited literature. This type of correlator has been in use for years
in the analysis of dynamic light
scattering :cite:`schatzel88a`. About a decade later it
found its way to the Fluorescence Correlation Spectroscopy
(FCS) :cite:`magatti01a`. The book of Frenkel and
Smit :cite:`frenkel02b` describes its application for the
special case of the velocity autocorrelation function.

.. _fig_correlator_scheme:

.. figure:: figures/correlator_scheme.png
   :scale: 50 %
   :alt: Schematic representation of buffers in the correlator.

   Schematic representation of buffers in the correlator.

Let us consider a set of :math:`N` observable values as schematically
shown in the figure above, where a value of index :math:`i` was
measured at times :math:`i\delta t`. We are interested in computing the
correlation function for a range 
of lag times :math:`\tau = (i-j)\delta t` between the measurements 
:math:`i` and :math:`j`. To simplify the notation, we drop
:math:`\delta t` when referring to observables and lag times.

The trivial implementation takes all possible pairs of values
corresponding to lag times
:math:`\tau \in [{\tau_{\mathrm{min}}}:{\tau_{\mathrm{max}}}]`. Without
loss of generality, we consider
:math:`{\tau_{\mathrm{min}}}=0`. The computational effort for such an
algorithm scales as
:math:`{\cal O} \bigl({\tau_{\mathrm{max}}}^2\bigr)`. As a rule of
thumb, this is feasible if :math:`{\tau_{\mathrm{max}}}< 10^3`. The
multiple tau correlator provides a solution to compute the correlation
functions for arbitrary range of the lag times by coarse-graining the
high :math:`\tau` values. It applies the naive algorithm to a relatively
small range of lag times :math:`\tau \in [0:p-1]` 
(:math:`p` corresponds to parameter ``tau_lin``). 
This we refer to as compression level 0. 
To compute the correlations for lag times
:math:`\tau \in [p:2(p-1)]`, the original data are first coarse-grained,
so that :math:`m` values of the original data are compressed to produce
a single data point in the higher compression level. Thus the lag time
between the neighboring values in the higher compression level
increases by a factor of :math:`m`, while the number of stored values
decreases by the same factor and the number of correlation operations at
this level reduces by a factor of :math:`m^2`. Correlations for lag
times :math:`\tau \in [2p:4(p-1)]` are computed at compression level 2,
which is created in an analogous manner from level 1. This can continue
hierarchically up to an arbitrary level for which enough data is
available. Due to the hierarchical reduction of the data, the algorithm
scales as
:math:`{\cal O} \bigl( p^2 \log({\tau_{\mathrm{max}}}) \bigr)`. Thus an
additional order of magnitude in :math:`{\tau_{\mathrm{max}}}` costs
just a constant extra effort.

The speedup is gained at the expense of statistical accuracy. The loss
of accuracy occurs at the compression step. In principle one can use any
value of :math:`m` and :math:`p` to tune the algorithm performance.
However, it turns out that using a high :math:`m` dilutes the data at
high :math:`\tau`. Therefore :math:`m=2` is hard-coded in the correlator
and cannot be modified by user. The value of :math:`p` remains an
adjustable parameter which can be modified by user by setting when
defining a correlation. In general, one should choose :math:`p \gg m` to
avoid loss of statistical accuracy. Choosing :math:`p=16` seems to be
safe but it may depend on the properties of the analyzed correlation
functions. A detailed analysis has been performed in
Ref. :cite:`ramirez10a`.

The choice of the compression function also influences the statistical
accuracy and can even lead to systematic errors. The default compression
function is which discards the second for the compressed values and
pushes the first one to the higher level. This is robust and can be
applied universally to any combination of observables and correlation
operation. On the other hand, it reduces the statistical accuracy as the
compression level increases. In many cases, the compression operation
can be applied, which averages the two neighboring values and the
average then enters the higher level, preserving almost the full
statistical accuracy of the original data. In general, if averaging can
be safely used or not, depends on the properties of the difference

.. math::

   \frac{1}{2} (A_i \otimes B_{i+p} + A_{i+1} \otimes B_{i+p+1} ) -
   \frac{1}{2} (A_i + A_{i+1} ) \otimes \frac{1}{2} (B_{i+p} +  B_{i+p+1})
   \label{eq:difference}

For example in the case of velocity autocorrelation function, the
above-mentioned difference has a small value and a random sign,  
different contributions cancel each other. On the other hand, in the of
the case of mean square displacement the difference is always positive,
resulting in a non-negligible systematic error. A more general
discussion is presented in Ref. :cite:`ramirez10a`.

.. _Accumulators:

Accumulators
------------

.. _Mean-variance calculator:

Mean-variance calculator
~~~~~~~~~~~~~~~~~~~~~~~~

In order to calculate the running mean and variance of an observable
:class:`espressomd.accumulators.MeanVarianceCalculator` can be used::

    import espressomd
    import espressomd.observables
    import espressomd.accumulators

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    system.part.add(id=0, pos=[5.0, 5.0, 5.0])
    position_observable = espressomd.observables.ParticlePositions(ids=(0,))
    accumulator = espressomd.accumulators.MeanVarianceCalculator(
        obs=position_observable, delta_N=1)
    system.auto_update_accumulators.add(accumulator)
    # Perform integration (not shown)
    print accumulator.get_mean()
    print accumulator.get_variance()

In the example above the automatic update of the accumulator is used. However,
it's also possible to manually update the accumulator by calling
:meth:`espressomd.accumulators.MeanVarianceCalculator.update`.

Cluster analysis
----------------

|es| provides support for online cluster analysis. Here, a cluster is a group of particles, such that you can get from any particle to any second particle by at least one path of neighboring particles.
I.e., if particle B is a neighbor of particle A, particle C is a neighbor of A and particle D is a neighbor of particle B, all four particles are part of the same cluster.
The cluster analysis is available in parallel simulations, but the analysis is carried out on the head node, only.


Whether or not two particles are neighbors is defined by a pair criterion. The available criteria can be found in :mod:`espressomd.pair_criteria`.
For example, a distance criterion which will consider particles as neighbors if they are closer than 0.11 is created as follows::

    from espressomd.pair_criteria import DistanceCriterion
    dc = DistanceCriterion(cut_off=0.11)

To obtain the cluster structure of a system, an instance of :class:`espressomd.cluster_analysis.ClusterStructure` has to be created.
To to create a cluster structure with above criterion::

    from espressomd.cluster_analysis import ClusterStructure
    cs = ClusterStructure(distance_criterion=dc)

In most cases, the cluster analysis is carried out by calling the :any:`espressomd.cluster_analysis.ClusterStructure.run_for_all_pairs` method. When the pair criterion is purely based on bonds,  :any:`espressomd.cluster_analysis.ClusterStructure.run_for_bonded_particles` can be used.

The results can be accessed via ClusterStructure.clusters, which is an instance of
:any:`espressomd.cluster_analysis.Clusters`.


Individual clusters are represented by instances of
:any:`espressomd.cluster_analysis.Cluster`, which provides access to the particles contained in a cluster as well as per-cluster analysis routines such as radius of gyration, center of mass and longest distance.
Note that the cluster objects do not contain copies of the particles, but refer to the particles in the simulation. Hence, the objects become outdated if the simulation system changes. On the other hand, it is possible to directly manipulate the particles contained in a cluster.





