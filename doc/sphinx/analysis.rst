.. _Analysis:

Analysis
========

The :mod:`espressomd.analyze` module provides online-calculation of local and global observables.
These exist for convenience and historical reasons.
Since arithmetic in TCL was quite painful to perform, there exists a series of common analysis routines which can be used whenever possible.
They usually have parts of the calculations performed in the core.


On the other hand, some observables are computed and stored in the
C-core of during a call to the function , while they are set up and
their results are collected from the script level. These observables are
more complex to implement and offer less flexibility, while they are
significantly faster and more memory efficient, and they can be set up
to be computed every few time steps. The observables in this class are
described in chapter :ref:`Analysis in the core`:
.

.. _Available observables:

Available observables
---------------------

.. _Minimal distances between particles:

Minimal distances between particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`espressomd.analyze.Analysis.mindist`
Returns the minimal distance between all particles in the system.

When used with type-lists as arguments, then the minimal distance between particles of only those types is determined.


:meth:`espressomd.analyze.Analysis.dist_to()`

Returns the minimal distance of all particles to either a particle (when used with an argument `id`) 
or a position coordinate when used with a vector `pos`.

For example, ::

    >>> import espressomd
    >>> system = espressomd.System()
    >>> system.box_l = [100, 100, 100]
    >>> for i in range(10):
    >>>     system.part.add(id=i, pos=[1.0, 1.0, i**2], type=0)
    >>> system.analysis.dist_to(id=4)
    7.0
    >>> system.analysis.dist_to(pos=[0,0,0])
    1.4142135623730951
    >>> system.analysis.mindist()
    1.0
    

.. _Particles in the neighbourhood:

Particles in the neighbourhood
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`espressomd.analyze.Analysis.nbhood`
 
Returns a list of the particle ids of that fall within a given radius of a target position.

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

    >>> system = espressomd.System()
    >>> box_l=10.
    >>> system.box_l = [box_l, box_l, box_l]
    >>> for i in range(5):
    >>>     system.part.add(id=i, pos=i*system.box_l, type=0)
    >>> bins, count=system.analysis.distribution(type_list_a=[0], type_list_b=[0], r_min=0.0, r_max = 10.0, r_bins=10)
    >>>
    >>> print(bins)
    [ 0.5  1.5  2.5  3.5  4.5  5.5  6.5  7.5  8.5  9.5]
    >>> print(count)
    [ 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]


.. _Radial density map:

Radial density map
~~~~~~~~~~~~~~~~~~
.. todo:: This feature is not implemented

analyze radial\_density\_map

Returns the radial density of particles around a given axis. Parameters
are:

-  histogram bins in x direction.

-  histogram bins in y direction.

-  range for analysis in x direction.

-  range for analysis in y direction.

-  rotate around given axis. (x, y, or z)

-  rotate around given point.

-  only analyze beads of given types.

-  histogram bins in angle theta.

This command does not do what you might expect. Here is an overview of
the currently identified properties.

#. is the number of bins along the axis of rotation.

#. is the number of bins in the radial direction.

#. The centre point () of the cylinder is located in the lower cap,
   i.e., is the height of the cylinder with respect to this centre
   point.

#. The bins are distributed along starting from 0 ().

#. The seem to average with respect to the centre of mass of the
   particles in the individual bins rather than with respect to the
   central axis, which one would think is natural.


.. _Cylindrical average:

Cylindrical Average
~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.cylindrical_average`

Calculates the particle distribution using cylindrical binning.

The volume considered is inside a cylinder defined by the parameters `center`, `axis`, `length` and  `radius`.

The geometrical details of the cylindrical binning is defined using ` bins_axial` and `bins_radial` which are the number bins in the axial and radial directions (respectively).
See figure :ref:`cylindrical_average` for a visual representation of the binning geometry.

.. _cylindrical_average:

.. figure:: figures/analysis_cylindrical_average.png
   :alt: Geometry for the cylindrical binning
   :align: center
   :height: 6.00000cm

   Geometry for the cylindrical binning


The command returns a list of lists. The outer list contains all data
combined whereas each inner list contains one line. Each lines stores a
different combination of the radial and axial index. The output might
look something like this

::

    [ [ 0 0 0.05 -0.25 0.0314159 0 0 0 0 0 0 ]
      [ 0 1 0.05 0.25 0.0314159 31.831 1.41421 1 0 0 0 ]
      ... ]

In this case two different particle types were present.
The columns of the respective lines are coded like this

=============    ============  ===========  ==========  =========  =======  ========   ========  =======  =========  =======
index_radial     index_axial   pos_radial   pos_axial   binvolume  density  v_radial   v_axial   density  v_radial   v_axial 
=============    ============  ===========  ==========  =========  =======  ========   ========  =======  =========  =======
0                0             0.05         -0.25       0.0314159  0        0          0         0        0          0      
0                1             0.05         0.25        0.0314159  31.831   1.41421    1         0        0          0      
=============    ============  ===========  ==========  =========  =======  ========   ========  =======  =========  =======

As one can see the columns `density`, `v_radial` and `v_axial` appear twice.
The order of appearance corresponds to the order of the types in the argument `types`.
For example if was set to `types=[0, 1]` then the first triple is associated to type 0 and
the second triple to type 1.

.. _Modes:

.. _Vkappa:

Vkappa
~~~~~~
:meth:`espressomd.analyze.Analysis.v_kappa`

.. todo:: Implementation appears to be incomplete

Calculates the compressibility :math:`V \times \kappa_T` through the
Volume fluctuations
:math:`V \times \kappa_T = \beta \left(\langle V^2\rangle - \langle V \rangle^2\right)`
:cite:`kolb99a`. Given no arguments this function calculates
and returns the current value of the running average for the volume
fluctuations.The `mode=reset` argument clears the currently stored values. With `mode=read` the
cumulative mean volume, cumulative mean squared volume and how many
samples were used can be retrieved. Likewise the option `mode=set` enables you to
set those.


.. _Radial distribution function:

Radial distribution function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.rdf`

Calculates a radial distribution function.


.. _Structure factor:

Structure factor
~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.structure_factor`

Calculate the structure factor for given types.

Returns the spherically averaged structure factor :math:`S(q)` of
particles specified in . :math:`S(q)` is calculated for all possible
wave vectors, :math:`\frac{2\pi}{L} <= q <= \frac{2\pi}{L}` `order`.


.. _Van-Hove autocorrelation function:

Van-Hove autocorrelation function :math:`G(r,t)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. todo:: This feature is not implemented

analyze vanhove

Returns the van Hove auto correlation function :math:`G(r,t)` and the
mean square displacement :math:`msd(t)` for particles of type for the
configurations stored in the array configs. This tool assumes that the
configurations stored with (see section ) are stored at equidistant time
intervals. :math:`G(r,t)` is calculated for each multiple of this time
intervals. For each time t the distribution of particle displacements is
calculated according to the specification given by , and . Optional
argument defines the maximum value of :math:`t` for which :math:`G(r,t)`
is calculated. If it is omitted or set to zero, maximum possible value
is used. If the particles perform a random walk (a normal diffusion
process) :math:`G(r,t)/r^2` is a Gaussian distribution for all times.
Deviations of this behavior hint on another diffusion process or on the
fact that your system has not reached the diffusive regime. In this case
it is also very questionable to calculate a diffusion constant from the
mean square displacement via the Stokes-Einstein relation.

The output corresponds to the blockfile format (see section ):

{ msd { …} } { vanhove { { …} { …} } }

The :math:`G(r,t)` are normalized such that the integral over space
always yields :math:`1`.


.. _Center of mass:

Center of mass
~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.center_of_mass`

Returns the center of mass of particles of the given type given by `part_type`.


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


.. _Aggregation:

Aggregation
~~~~~~~~~~~
.. todo:: This feature is not implemented

analyze aggregation

Returns the aggregate size distribution for the molecules in the
molecule id range to . If any monomers in two different molecules are
closer than they are considered to be in the same aggregate. One can use
the optional parameter to specify a minimum number of contacts such that
only molecules having at least contacts will be considered to be in the
same aggregate. The second optional parameter enables one to consider
aggregation state of only oppositely charged particles.


.. _Identifying pearl necklace structures:


.. _Finding holes:

Finding holes
~~~~~~~~~~~~~
.. todo:: This feature is not implemented

analyze holes

Function for the calculation of the unoccupied volume (often also called
free volume) in a system. Details can be found in
:cite:`schmitz00a`. It identifies free space in the
simulation box via a mesh based cluster algorithm. Free space is defined
via a probe particle and its interactions with other particles which
have to be defined through LJ interactions with the other existing
particle types via the inter command before calling this routine. A
point of the mesh is counted as free space if the distance of the point
is larger than LJ\_cut+LJ\_offset to any particle as defined by the LJ
interaction parameters between the probe particle type and other
particle types.How to use this function: Define interactions between all
(or the ones you are interested in) particle types in your system and a
fictitious particle type. Practically one uses the van der Waals radius
of the particles plus the size of the probe you want to use as the
Lennard Jones cutoff. The mesh spacing is the box length divided by the
.

{ { } { } { } }

A hole is defined as a continuous cluster of mesh elements that belong
to the unoccupied volume. Since the function is quite rudimentary it
gives back the whole information suitable for further processing on the
script level. and are given in number of mesh points, which means you
have to calculate the actual size via the corresponding volume or
surface elements yourself. The complete information is given in the
element\_lists for each hole. The element numbers give the position of a
mesh point in the linear representation of the 3D grid (coordinates are
in the order x, y, z). Attention: the algorithm assumes a cubic box.
Surface results have not been tested. .


.. _Temperature of the lb fluid:

Temperature of the LB fluid
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. todo:: This feature is not implemented

This command returns the temperature of the lattice-Boltzmann (LB)
fluid, see Chapter [sec:lb], by averaging over the fluid nodes. In case
or are compiled in and boundaries are defined, only the available fluid
volume is taken into account.


.. _Momentum of the system:

Momentum of the System
~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.analyze_linear_momentum`

This command returns the total linear momentum of the particles and the
lattice-Boltzmann (LB) fluid, if one exists. Giving the optional
parameters either causes the command to ignore the contribution of LB or
of the particles.


.. _Energies:

Energies
~~~~~~~~
:meth:`espressomd.analyze.Analysis.energy`

Returns the energies of the system.
The the different energetic contributions to the total energy can also be obtained (kinetic, bonded,non-bonded, coublomb)).


.. _Pressure:

Pressure
~~~~~~~~

:meth:`espressomd.analyze.Analysis.pressure`

Computes the pressure and its contributions in the system. It
returns all the contributions to the total pressure (see :meth:`espressomd.analyze.Analysis.pressure`).

The pressure is calculated (if there are no electrostatic interactions)
by

.. math::
     p = \frac{2E_{kinetic}}{Vf} + \frac{\sum_{j>i} {F_{ij}r_{ij}}}{3V}
     :label: eqptens

where :math:`f=3` is the number of translational degrees of freedom of
each particle, :math:`V` is the volume of the system,
:math:`E_{kinetic}` is the kinetic energy, :math:`F_{ij}` the force
between particles i and j, and :math:`r_{ij}` is the distance between
them. The kinetic energy divided by the degrees of freedom is

.. math:: \frac{2E_{kinetic}}{f} = \frac{1}{3}\sum_{i} {m_{i}v_{i}^{2}}.

Note that Equation :eq:`eqptens` can only be applied to pair potentials and
central forces. Description of how contributions from other interactions
are calculated is beyond the scope of this manual. Three body potentials
are implemented following the procedure in
Ref. :cite:`thompson09a`. A different formula is used to
calculate contribution from electrostatic interactions in P3M. For
electrostatic interactions, the :math:`k`-space contribution is not well
tested, so use with caution! Anything outside that is currently not
implemented. Four-body dihedral potentials are not included. Except of 
VIRTUAL\_SITES\_RELATIVE constraints all other
constraints of any kind are not currently accounted for in the pressure
calculations. The pressure is no longer correct, e.g., when particles
are confined to a plane.

The command is implemented in parallel.

.. _Stress Tensor:

Stress Tensor
~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.stress_tensor`

Computes the stress tensor of the system with options which are
described by in :meth: espressomd.System.analysis.stress_tensor. 
It is called a stress tensor but the sign convention follows that of a pressure tensor.

The virial stress tensor is calculated by

.. math:: p^{(kl)} = \frac{\sum_{i} {m_{i}v_{i}^{(k)}v_{i}^{(l)}}}{V} + \frac{\sum_{j>i}{F_{ij}^{(k)}r_{ij}^{(l)}}}{V}

where the notation is the same as for in and the superscripts :math:`k`
and :math:`l` correspond to the components in the tensors and vectors.

Note that the angular velocities of the particles are not included in
the calculation of the stress tensor.

The command is implemented in parallel.


.. _Local Stress Tensor:

Local Stress Tensor
~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.analyze.Analysis.local_stress_tensor`
.. todo:: This feature is not implemented

A cuboid is defined in the system and divided into bins.
For each of these bins a stress tensor is calculated using the Irving Kirkwood method.
That is, a given interaction contributes towards the stress tensor in a bin proportional to the fraction of the line connecting the two particles that is within the bin.

If the P3M and MMM1D electrostatic methods are used, these interactions
are not included in the local stress tensor. The DH and RF methods, in
contrast, are included. Concerning bonded interactions only two body
interactions (FENE, Harmonic) are included (angular and dihedral are
not). For all electrostatic interactions only the real space part is
included.

Care should be taken when using constraints of any kind, since these are
not accounted for in the local stress tensor calculations.

The command is implemented in parallel.

{ { LocalStressTensor } { { } { } } }

specifying the local pressure tensor in each bin.



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
vector of its centre of mass. The sum runs over all :math:`N` particles
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


.. _Internal distances:

Internal distances
^^^^^^^^^^^^^^^^^^
.. todo:: This feature is not implemented

analyze

Returns the averaged internal distances within the chains (over all
pairs of particles). If is used, the values are averaged over all stored
configurations (see section ).

{ … }

The index corresponds to the number of beads between the two monomers
considered (0 = next neighbours, 1 = one monomer in between, …).


.. _Internal distances II (specific monomer):

Internal distances II (specific monomer)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. todo:: This feature is not implemented

analyze

In contrast to , it does not average over the whole chain, but rather
takes the chain monomer at position (default: :math:`0`, the first
monomer on the chain) to be the reference point to which all internal
distances are calculated. If is used, the values will be averaged over
all stored configurations (see section ).

{ … }


.. _Bond lengths:

Bond lengths
^^^^^^^^^^^^
.. todo:: This feature is not implemented

analyze

Analyzes the bond lengths of the chains in the system. Returns its
average, the standard deviation, the maximum and the minimum. If you
want to look only at specific chains, use the optional arguments,
:math:`\var{chain\_start} =
2*\var{MPC}` and :math:`\var{n\_chains} = 1` to only include the third
chain’s monomers. If is used, the value will be averaged over all stored
configurations (see section ). This function assumes linear chain
topology and does not check if the bonds really exist!

{ }


.. _Form factor:

Form factor
^^^^^^^^^^^
.. todo:: This feature is not implemented

| analyze

Computes the spherically averaged form factor of a single chain, which
is defined by

.. math::

   S(q) = \frac{1}{\var{chain\_length}} \sum_{i,j=1}^{\var{chain\_length}}
     \frac{\sin(q r_{ij})}{q r_{ij}}

of a single chain, averaged over all chains for :math:`\var{qbin}+1`
logarithmically spaced q-vectors :math:`\var{qmin}, \dots ,\var{qmax}`
where :math:`\var{qmin}>0` and :math:`\var{qmax}>\var{qmin}`. If is
used, the form factor will be averaged over all stored configurations
(see section ).

{ { } }

with :math:`q \in \{\var{qmin},\dots,\var{qmax}\}`.


.. _Chain radial distribution function:

Chain radial distribution function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:meth:`espressomd.analyze.Analysis.rdf_chain`

Returns three radial distribution functions (rdf) for the chains.
The first rdf is calculated for monomers belonging to different chains,
the second rdf is for the centers of mass of the chains and 
the third one is the distribution of the closest distances between the chains (the
shortest monomer-monomer distances).


