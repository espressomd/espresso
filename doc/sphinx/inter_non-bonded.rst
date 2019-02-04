.. _Non-bonded interactions:

Non-bonded interactions
=======================

In |es|, interactions are set up and investigated by the :mod:`espressomd.interactions` module. There are
mainly two types of interactions: non-bonded and bonded interactions.

Non-bonded interactions only depend on the *type* of the two particles
involved. This also applies to the electrostatic interaction; however,
due to its long-ranged nature, it requires special care and |es| handles it
separately with a number of state-of-the-art algorithms. To specify particle
type and charge see :ref:`Setting up particles`.

A bonded interaction defines an interaction between a number of specific
particles; it only applies to the set of particles for which it has been
explicitly set. A bonded interaction between a set of particles has to
be specified explicitly by the command, while the command is used to
define the interaction parameters.

.. todo::
    IMPLEMENT: print interaction list

.. _Isotropic non-bonded interactions:

Isotropic non-bonded interactions
---------------------------------

Non-bonded interaction are configured via the :class:`espressomd.interactions.NonBondedInteraction` class, which is a member of :class:`espressomd.system.System`::

    system.non_bonded_inter[type1, type2]

This command defines an interaction between all particles of type ``type1`` and
``type2``. Possible interaction types and their parameters are
listed below.

.. todo::
    Implement this functionality:
    If the interaction is omitted, the command returns the
    currently defined interaction between the two types using the syntax to
    define the interaction

For many non-bonded interactions, it is possible to artificially cap the
forces, which often allows to equilibrate the system much faster. See
the subsection :ref:`Capping the force during warmup` for more details.

.. _Tabulated interaction:

Tabulated interaction
~~~~~~~~~~~~~~~~~~~~~

.. note ::

    Feature ``TABULATED`` required.


The interface for tabulated interactions are implemented in the
:class:`~espressomd.interactions.TabulatedNonBonded` class. They can be configured
via the following syntax::

  system.non_bonded_inter[type1, type2].tabulated.set_params(
      min='min', max='max', energy='energy', force='force')


This defines an interaction between particles of the types ``type1`` and
``type2`` according to an arbitrary tabulated pair potential by linear interpolation.
``force`` specifies the tabulated forces and ``energy`` the energies as a function of the
separation distance. ``force`` and ``energy`` have to have the same length :math:`N_\mathrm{points}`.
Take care when choosing the number of points, since a copy of each lookup
table is kept on each node and must be referenced very frequently.
The maximal tabulated separation distance also acts as the effective cutoff
value for the potential.

The values of :math:`r` are assumed to be equally distributed between
:math:`r_\mathrm{min}` and :math:`r_\mathrm{max}` with a fixed distance
of :math:`(r_\mathrm{max}-r_\mathrm{min})/(N_\mathrm{points}-1)`.

.. _Lennard-Jones interaction:

Lennard-Jones interaction
~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``LENNARD_JONES`` required.

The interface for the Lennard-Jones interaction is implemented in
:class:`~espressomd.interactions.LennardJonesInteraction`. The Lennard-Jones parameters
can be set via::

    system.non_bonded_inter[type1, type2].lennard_jones.set_params(**kwargs)

This command defines the traditional (12-6)-Lennard-Jones interaction
between particles of the types ``type1`` and ``type2``. For a description of the input arguments
see :class:`~espressomd.interactions.LennardJonesInteraction`. The potential is defined by

.. math::

   \label{eq:lj}
     V_\mathrm{LJ}(r) =
       \begin{cases}
         4 \epsilon \left[ \left(\frac{\sigma}{r-r_\mathrm{off}}\right)^{12}
         - \left(\frac{\sigma}{r-r_\mathrm{off}}\right)^6+c_\mathrm{shift}\right]
         & \mathrm{if~} r_\mathrm{min}+r_\mathrm{off} < r < r_\mathrm{cut}+r_\mathrm{off}\\
         0
         & \mathrm{otherwise}
       \end{cases}.

The traditional Lennard-Jones potential is the "work-horse" potential of
particle--particle interactions in coarse-grained simulations. It is a
simple model for the van-der-Waals interaction, and is attractive at
large distance, but strongly repulsive at short distances.
:math:`r_\mathrm{off} + \sigma` corresponds to the sum of
the radii of the interaction particles. At this distance, the potential is
:math:`V_\mathrm{LJ}(r_\mathrm{off} + \sigma) = 4 \epsilon c_\mathrm{shift}`.
The minimum of the potential is at
:math:`V_\mathrm{LJ}(r_\mathrm{off} +
2^\frac{1}{6}\sigma) =
-\epsilon + 4 \epsilon c_\mathrm{shift}`. Beyond this value the interaction is attractive.
Beyond the distance :math:`r_\mathrm{cut}` the potential is cut off and the interaction force is zero.

If :math:`c_\mathrm{shift}` is not set or it is set to the string *auto*, the shift will be
automatically computed such that the potential is continuous at the
cutoff radius. If is not set, it is set to :math:`0`.

The net force on a particle can be capped by using force capping , see
section :ref:`Capping the force during warmup`

An optional additional parameter can be used to restrict the interaction
from a *minimal* distance :math:`r_\mathrm{min}`. This is an
optional parameter, set to :math:`0` by default.

A special case of the Lennard-Jones potential is the
Weeks-Chandler-Andersen (WCA) potential, which one obtains by putting
the cutoff into the minimum, choosing
:math:`r_\mathrm{cut}=2^\frac{1}{6}\sigma`. The WCA
potential is purely repulsive, and is often used to mimic hard sphere
repulsion.

.. _Generic Lennard-Jones interaction:

Generic Lennard-Jones interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``LENNARD_JONES_GENERIC`` required.


The interface for the generic Lennard-Jones interactions is implemented in
:class:`espressomd.interactions.GenericLennardJonesInteraction`. They
are configured via the syntax::

    system.non_bonded_inter[type1, type2].generic_lennard_jones.set_params(**kwargs)

This command defines a generalized version of the Lennard-Jones
interaction (see :ref:`Lennard-Jones interaction`) between particles of the
types ``type1`` and ``type2``. The potential is defined by

.. math::

   \label{eq:lj-generic}
     V_\mathrm{LJ}(r) =
       \begin{cases}
         \epsilon\left[b_1\left(\frac{\sigma}{r-r_\mathrm{off}}\right)^{e_1}
         -b_2\left(\frac{\sigma}{r-r_\mathrm{off}}\right)^{e_2}+c_\mathrm{shift}\right]
         & \mathrm{if~} r_\mathrm{min}+r_\mathrm{off} < r < r_\mathrm{cut}+r_\mathrm{off}\\
         0
         & \mathrm{otherwise}
       \end{cases}\ .

Note that the prefactor 4 of the standard LJ potential is missing, so
the normal LJ potential is recovered for :math:`b_1=b_2=4`,
:math:`e_1=12` and :math:`e_2=6`.

The net force on a particle can be capped by using force capping ``system.non_bonded_inter.set_force_cap(max)``, see
section :ref:`Capping the force during warmup`

The optional ``LJGEN_SOFTCORE`` feature activates a softcore version of
the potential, where the following transformations apply:
:math:`\epsilon \rightarrow \lambda \epsilon` and
:math:`r-r_\mathrm{off} \rightarrow \sqrt{(r-r_\mathrm{off})^2 +
(1-\lambda) \delta \sigma^2}`. :math:`\lambda` allows to tune the strength of the
interaction, while :math:`\delta` varies how smoothly the potential goes to zero as
:math:`\lambda\rightarrow 0`. Such a feature allows one to perform
alchemical transformations, where a group of atoms can be slowly turned
on/off during a simulation.

.. _Weeks-Chandler-Andersen interaction:

Weeks-Chandler-Andersen interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``WCA`` required.


The interface for the Weeks-Chandler-Andersen interactions is implemented in
:class:`espressomd.interactions.WCAInteraction`. They
are configured via the syntax::

    system.non_bonded_inter[type1, type2].wca.set_params(**kwargs)

This command defines a Weeks-Chandler-Andersen interaction between particles of the
types ``type1`` and ``type2``. The potential is defined by

.. math::

   \label{eq:wca}
     V_\mathrm{WCA}(r) =
       \begin{cases}
         4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12}
         - \left(\frac{\sigma}{r}\right)^6 + \frac{1}{4} \right]
         & \mathrm{if~} r < \sigma 2^{\frac{1}{6}}\\
         0
         & \mathrm{otherwise}
       \end{cases}.

The net force on a particle can be capped by using
force capping ``system.non_bonded_inter.set_force_cap(max)``, see
section :ref:`Capping the force during warmup`

.. _Lennard-Jones cosine interaction:

Lennard-Jones cosine interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Feature ``LJCOS`` and/or ``LJCOS2`` required.

.. code::

   system.non_bonded_inter[type1, type2].lennard_jones_cos.set_params(**kwargs)
   system.non_bonded_inter[type1, type2].lennard_jones_cos2.set_params(**kwargs)

:class:`espressomd.interactions.LennardJonesCosInteraction` and
:class:`espressomd.interactions.LennardJonesCos2Interaction` specifies
a Lennard-Jones interaction with cosine tail :cite:`soddeman01a`
between particles of the types ``type1`` and ``type2``. The first variant
behaves as follows: Until the minimum of the Lennard-Jones potential
at :math:`r_\mathrm{min} = r_\mathrm{off} + 2^{\frac{1}{6}}\sigma`, it
behaves identical to the unshifted Lennard-Jones potential
(:math:`c_\mathrm{shift}=0`). Between :math:`r_\mathrm{min}` and :math:`r_\mathrm{cut}`, a cosine is used to
smoothly connect the potential to 0, i.e.,

.. math::

    V(r)=\frac{1}{2}\epsilon\left(\cos\left[\alpha(r - r_\mathrm{off})^2 + \beta\right]-1\right),

where :math:`\alpha = \pi\left[(r_\mathrm{cut} -
r_\mathrm{off})^2-(r_\mathrm{min} - r_\mathrm{off})^2\right]^{-1}` and
:math:`\beta = \pi - \left(r_\mathrm{min} -
r_\mathrm{off}\right)^2\alpha`.

In the second variant, the cutoff radius is
:math:`r_\mathrm{cut}=r_\mathrm{min} + \omega`, where
:math:`r_\mathrm{min} = r_\mathrm{off} + 2^{\frac{1}{6}}\sigma` as in
the first variant. The potential between :math:`r_\mathrm{min}` and
:math:`r_\mathrm{cut}` is given by

.. math::

   V(r)=-\epsilon\cos^2\left[\frac{\pi}{2\omega}(r - r_\mathrm{min})\right].

For :math:`r < r_\mathrm{min}`, :math:`V(r)` is implemented as normal
:ref:`Lennard-Jones interaction` with :math:`c_\mathrm{shift} = 0`.

The net force on a particle can be capped by using force capping, see
section :ref:`Capping the force during warmup`

.. _Smooth step interaction:

Smooth step interaction
~~~~~~~~~~~~~~~~~~~~~~~

.. note::
     Feature ``SMOOTH_STEP`` required.

The interface for the smooth-step interaction is implemented in
:class:`espressomd.interactions.SmoothStepInteraction`. The smooth-step parameters
can be set via::

     system.non_bonded_inter[type1, type2].smooth_step.set_params(**kwargs)

This defines a smooth step interaction between particles of the types ``type1``
and ``type2``, for which the potential is

.. math:: V(r)= \left(d/r\right)^n + \epsilon/(1 + \exp\left[2k_0 (r - \sigma)\right])

for :math:`r<r_\mathrm{cut}`, and :math:`V(r)=0` elsewhere. With
:math:`n` around 10, the first term creates a short range repulsion
similar to the Lennard-Jones potential, while the second term provides a
much softer repulsion. This potential therefore introduces two length
scales, the range of the first term, :math:`d`, and the range of
the second one, :math:`\sigma`, where in general :math:`d<\sigma`.

.. _BMHTF potential:

BMHTF potential
~~~~~~~~~~~~~~~

.. note::
     Feature ``BMHTF_NACL`` required.

The interface for the smooth-step interaction is implemented in
:class:`espressomd.interactions.BMHTFInteraction`. The parameters of the BMHTF potential
can be set via::

     system.non_bonded_inter[type1, type2].bmhtf.set_params(**kwargs)

This defines an interaction with the *short-ranged part* of the
Born-Meyer-Huggins-Tosi-Fumi potential between particles of the types ``type1``
and ``type2``, which is often used to simulate NaCl crystals. The potential is
defined by:

.. math::

   V(r)= A\exp\left[B(\sigma - r)\right] -
     C r^{-6} - D r^{-8} + \epsilon_\mathrm{shift},

where :math:`\epsilon_\mathrm{shift}` is automatically chosen such that
:math:`V(r_\mathrm{cut})=0`. For
:math:`r\ge r_\mathrm{cut}`, the :math:`V(r)=0`.

For NaCl, the parameters should be chosen as follows:

+---------+---------------------------------------------------------+-----------------------------------------------------------+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------------------------------------+
| types   | :math:`A` :math:`\left(\mathrm{kJ}/\mathrm{mol}\right)` | :math:`B` :math:`\left(\mathring{\mathrm{A}}^{-1}\right)` | :math:`C` :math:`\left(\mathring{\mathrm{A}}^6 \mathrm{kJ}/\mathrm{mol})\right)` | :math:`D` :math:`\left(\mathring{\mathrm{A}}^8 \mathrm{kJ}/\mathrm{mol}\right)` | :math:`\sigma` :math:`\left(\mathring{\mathrm{A}}\right)` |
+=========+=========================================================+===========================================================+==================================================================================+=================================================================================+===========================================================+
| Na-Na   | 25.4435                                                 | 3.1546                                                    | 101.1719                                                                         | 48.1771                                                                         | 2.34                                                      |
+---------+---------------------------------------------------------+-----------------------------------------------------------+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------------------------------------+
| Na-Cl   | 20.3548                                                 | 3.1546                                                    | 674.4793                                                                         | 837.0770                                                                        | 2.755                                                     |
+---------+---------------------------------------------------------+-----------------------------------------------------------+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------------------------------------+
| Cl-Cl   | 15.2661                                                 | 3.1546                                                    | 6985.6786                                                                        | 14031.5785                                                                      | 3.170                                                     |
+---------+---------------------------------------------------------+-----------------------------------------------------------+----------------------------------------------------------------------------------+---------------------------------------------------------------------------------+-----------------------------------------------------------+

The cutoff can be chosen relatively freely because the potential decays
fast; a value around 10 seems reasonable.

In addition to this short ranged interaction, one needs to add a
Coulombic, long-ranged part. If one uses elementary charges, a charge of
:math:`q=+1` for the Na-particles, and :math:`q=-1` for the
Cl-particles, the corresponding prefactor of the Coulomb interaction is
:math:`\approx 1389.3549\,\mathrm{kJ}/\mathrm{mol}`.

.. _Morse interaction:

Morse interaction
~~~~~~~~~~~~~~~~~

.. note::
     Feature ``MORSE`` required.

The interface for the Morse interaction is implemented in
:class:`espressomd.interactions.MorseInteraction`. The Morse interaction parameters
can be set via::

     system.non_bonded_inter[type1, type2].morse.set_params(**kwargs)

This defines an interaction using the Morse potential between particles
of the types ``type1`` and ``type2``. It serves similar purposes as the Lennard-Jones
potential, but has a deeper minimum, around which it is harmonic. This
models the potential energy in a diatomic molecule.

For :math:`r < r_\mathrm{cut}`, this potential is given by

.. math::

   V(r)=\epsilon\left(\exp\left[-2 \alpha \left(r - r_\mathrm{min}\right)\right]
       - 2\exp\left[-\alpha\left(r - r_\mathrm{min}\right)\right]\right) -
     \epsilon_\mathrm{shift},

where is again chosen such that :math:`V(r_\mathrm{cut})=0`. For
:math:`r\ge r_\mathrm{cut}`, the :math:`V(r)=0`.

.. _Buckingham interaction:

Buckingham interaction
~~~~~~~~~~~~~~~~~~~~~~

.. note::
     Feature ``BUCKINGHAM`` required.

The interface for the Buckingham interaction is implemented in
:class:`espressomd.interactions.BuckinghamInteraction`. The Buckingham interaction parameters
can be set via::

     system.non_bonded_inter[type1, type2].morse.set_params(**kwargs)

This defines a Buckingham interaction between particles of the types *type1* and *type2*,
for which the potential is given by

.. math:: V(r)= A \exp(-B r) - C r^{-6} - D r^{-4} + \epsilon_\mathrm{shift}

for :math:`r_\mathrm{discont} < r < r_\mathrm{cut}`. Below :math:`r_\mathrm{discont}`,
the potential is linearly continued towards :math:`r=0`, similarly to
force capping, see below. Above :math:`r=r_\mathrm{cut}`, the
potential is :math:`0`.

.. _Soft-sphere interaction:

Soft-sphere interaction
~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``SOFT_SPHERE`` required.

The interface for the Soft-sphere interaction is implemented in
:class:`espressomd.interactions.SoftSphereInteraction`. The Soft-sphere parameters
can be set via::

    system.non_bonded_inter[type1, type2].soft_sphere.set_params(**kwargs)

This defines a soft sphere interaction between particles of the types ``type1``
and ``type2``, which is defined by a single power law:

.. math:: V(r)=a\left(r-r_\mathrm{offset}\right)^{-n}

for :math:`r<r_\mathrm{cut}`, and :math:`V(r)=0` above. There is
no shift implemented currently, which means that the potential is
discontinuous at :math:`r=r_\mathrm{cut}`. Therefore energy
calculations should be used with great caution.

.. _Membrane-collision interaction:

Membrane-collision interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
     Feature ``MEMBRANE_COLLISION`` required.

This defines a membrane collision interaction between particles of the
types ``type1`` and ``type2``, where particle of ``type1`` belongs to one OIF or OIF-like object and
particle of ``type2`` belongs to another such object.

It is very similar to soft-sphere interaction, but it takes into account
the local outward normal vectors on the surfaces of the two objects to
determine the direction for repulsion of objects (i.e. determine whether
the two membranes are intersected). It is inversely proportional to the
distance of nodes of membranes that are not crossed and saturating with
growing distance of nodes of crossed membranes.

In order to work with the OIF objects, both OIF objects need to be created
using OifCellType class with keyword ``normal=1``, because this implicitly sets up the
bonded out-direction interaction, which computes the outward normal
vector.

The membrane-collision interaction for non-intersected membranes is then
defined by:

.. math:: V(d)= a\frac{1}{1+e^{n\left(d-d_\mathrm{offset}\right)}},

for :math:`d<d_\mathrm{cut}` and :math:`V(d)=0` above. For
intersected membranes, it is defined as :math:`V(-d)`. There is no shift
implemented currently, which means that the potential is discontinuous
at :math:`d=d_\mathrm{cut}`. Therefore energy calculations should
be used with great caution.

.. _Hat interaction:

Hat interaction
~~~~~~~~~~~~~~~

.. note::
    Feature ``HAT`` required.

The interface for the Lennard-Jones interaction is implemented in
:class:`espressomd.interactions.HatInteraction`. The hat parameters
can be set via::

    system.non_bonded_inter[type1, type2].hat.set_params(**kwargs)

This defines a simple force ramp between particles of two types.
The maximal force acts at zero distance and zero force is applied at
distances :math:`r_c` and bigger. For distances smaller than :math:`r_c`,
the force is given by

.. math:: F(r)=F_{\text{max}} \cdot \left( 1 - \frac{r}{r_c} \right),

for distances exceeding :math:`r_c`, the force is zero.

The potential energy is given by

.. math:: V(r)=F_{\text{max}} \cdot (r-r_c) \cdot \left( \frac{r+r_c}{2r_c} - 1 \right),

which is zero for distances bigger than :math:`r_c` and continuous at distance :math:`0`.

This is the standard conservative DPD potential and can be used in
combination with :ref:`Dissipative Particle Dynamics (DPD)`.



Hertzian interaction
~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``HERTZIAN`` required.

The interface for the Hertzian interaction is implemented in
:class:`espressomd.interactions.HertzianInteraction`. The Hertzian interaction parameters
can be set via::

    system.non_bonded_inter[type1, type2].hertzian.set_params(**kwargs)

This defines an interaction according to the Hertzian potential between
particles of the types ``type1`` and ``type2``. The Hertzian potential is defined by

.. math::

   V(r)=
     \begin{cases} \epsilon\left(1-\frac{r}{\sigma}\right)^{5/2} & r < \sigma\\
       0 & r \ge \sigma.
     \end{cases}

The potential has no singularity and is defined everywhere; the
potential has a non-differentiable maximum at :math:`r=0`, where the force
is undefined.

.. _Gaussian:

Gaussian
~~~~~~~~

.. note::
    Feature ``GAUSSIAN`` required.

The interface for the Gaussian interaction is implemented in
:class:`espressomd.interactions.GaussianInteraction`. The Gaussian interaction parameters
can be set via::

    system.non_bonded_inter[type1, type2].gaussian.set_params(**kwargs)

This defines an interaction according to the Gaussian potential between
particles of the types ``type1`` and ``type2``. The Gaussian potential is defined by

.. math::

   V(r) =
     \begin{cases} \epsilon\,e^{-\frac{1}{2}\left(\frac{r}{\sigma}\right)^{2}}
       & r < r_\mathrm{cut}\\
     0 & r \ge r_\mathrm{cut}
     \end{cases}

The Gaussian potential is smooth except at the cutoff, and has a finite
overlap energy of :math:`\epsilon`. It can be used to model overlapping
polymer coils.

Currently, there is no shift implemented, which means that the potential
is discontinuous at :math:`r=r_\mathrm{cut}`. Therefore use
caution when performing energy calculations. However, you can often
choose the cutoff such that the energy difference at the cutoff is less
than a desired accuracy, since the potential decays very rapidly.

.. _DPD interaction:

DPD interaction
~~~~~~~~~~~~~~~

.. note::
    Feature ``DPD`` required.

This is a special interaction that is to be used in conjunction with the
`Dissipative Particle Dynamics (DPD)` thermostat, for a general description
of the algorithm see there. The parameters can be set via::

    system.non_bonded_inter[type1, type2].dpd.set_params(**kwargs)

This command defines a velocity-dependent interaction
between particles of the types ``type1`` and ``type2``. For a description of the input arguments
see :class:`espressomd.interactions.DPDInteraction`. The interaction
only has an effect if the DPD thermostat activated, and acts according to the
temperature of the thermostat.

The interaction consists of a dissipative force :math:`\vec{F}_{ij}^{D}` and
a random force :math:`\vec{F}_{ij}^R`, and is decomposed into a component
parallel and perpendicular to the distance vector of the particle pair :math:`\vec{F}_{ij}`.
The parameters for the two parts can be chosen independently.
The force contributions of the parallel part are

.. math:: \vec{F}_{ij}^{D} = -\zeta w^D (r_{ij}) (\hat{r}_{ij} \cdot \vec{v}_{ij}) \hat{r}_{ij}

for the dissipative force and

.. math:: \vec{F}_{ij}^R = \sigma w^R (r_{ij}) \Theta_{ij} \hat{r}_{ij}

for the random force. Here :math:`w^D` and :math:`w^R` are weight functions that
can be specified via the weight_function parameter of the interaction. The dissipative
and random weight function are related by the dissipation-fluctuation theorem:

.. math:: (\sigma w^R (r_{ij}))^2=\zeta w^D (r_{ij}) \text{k}_\text{B} T

The possible values for weight_function are 0 and 1, corresponding to the
order of :math:`w^D`:

.. math::

   w^D (r_{ij}) = ( w^R (r_{ij})) ^2 =
      \left\{
   \begin{array}{clcr}
                1                      & , \; \text{weight_function} = 0 \\
                {( 1 - \frac{r_{ij}}{r_c}} )^2 & , \; \text{weight_function} = 1
      \end{array}
      \right.

For the perpendicular part, the dissipative force is calculated by

.. math:: \vec{F}_{ij}^{D} = -\zeta w^D (r_{ij}) (I-\hat{r}_{ij}\otimes\hat{r}_{ij}) \cdot \vec{v}_{ij}

The random force by

.. math:: \vec{F}_{ij}^R = \sigma w^R (r_{ij}) (I-\hat{r}_{ij}\otimes\hat{r}_{ij}) \cdot \vec{\Theta}_{ij}

The parameters define the strength of the friction and the cutoff in the
same way as above. Note: This interaction does *not* conserve angular
momentum.

.. _Thole correction:

Thole correction
~~~~~~~~~~~~~~~~

.. note::

    Requires features ``THOLE`` and ``ELECTROSTATICS``.

.. note::

    ``THOLE`` is only implemented for the P3M electrostatics solver.

The Thole correction is closely related to simulations involving
:ref:`Particle polarizability with thermalized cold Drude oscillators`.
In this context, it is used to correct for overestimation of
induced dipoles at short distances. Ultimately, it alters the short-range
electrostatics of P3M to result in a damped Coulomb interaction potential
:math:`V(r) = \frac{q_1 q_2}{r} \cdot (1- e^{-s r} (1 + \frac{s r}{2}) )`.  The
Thole scaling coefficient :math:`s` is related to the polarizabilities
:math:`\alpha` and Thole damping parameters :math:`a` of the interacting
species via :math:`s = \frac{ (a_i + a_j) / 2 }{ (\alpha_i \alpha_j)^{1/6} }`.
Note that for the Drude oscillators, the Thole correction should be applied
only for the dipole part :math:`\pm q_d` added by the Drude charge and not on
the total core charge, which can be different for polarizable ions. Also note
that the Thole correction acts between all dipoles, intra- and intermolecular.
Again, the accuracy is related to the P3M accuracy and the split between
short-range and long-range electrostatics interaction. It is configured by::

    system = espressomd.System()
    system.non_bonded_inter[type_1,type_2].thole.set_params(scaling_coeff=<float>, q1q2=<float>)

with parameters:
    * ``scaling_coeff``: The scaling coefficient :math:`s`.
    * ``q1q2``: The charge factor of the involved charges.

Because the scaling coefficient depends on the *mixed* polarizabilities and the
nonbonded interaction is controlled by particle types, each Drude charge with a
unique polarizability has to have a unique type. Each Drude charge type has
a Thole correction interaction with all other Drude charges and all Drude
cores, except the one it's connected to.  This exception is handled internally
by disabling Thole interaction between particles connected via Drude bonds.
Also, each Drude core has a Thole correction interaction with all other Drude
cores and Drude charges. To assist with the bookkeeping of mixed scaling
coefficients, the helper method :meth:`~espressomd.drude_helpers.add_drude_particle_to_core` (see
:ref:`Particle polarizability with thermalized cold Drude oscillators`)
collects all core types, Drude types and relevant parameters when a Drude
particle is created. The user already provided all the information when
setting up the Drude particles, so the simple call::

    add_all_thole(<system>, <verbose>)

given the :class:`espressomd.System() <espressomd.system.System>` object, uses this information to create all
necessary Thole interactions. The method calculates the mixed scaling
coefficient `s` and creates the non-bonded Thole interactions between the
collected types to cover all the Drude-Drude, Drude-core and core-core
combinations. No further calls of :meth:`~espressomd.drude_helpers.add_drude_particle_to_core` should
follow. Set ``verbose`` to ``True`` to print out the coefficients, charge factors
and involved types.

The samples folder contains the script :file:`drude_bmimpf6.py` with a fully
polarizable, coarse grained ionic liquid where this approach is applied.
To use the script, compile espresso with the following features:

.. code-block:: c++

    #define EXTERNAL_FORCES
    #define MASS
    #define LANGEVIN_PER_PARTICLE
    #define ROTATION
    #define ROTATIONAL_INERTIA
    #define ELECTROSTATICS
    #define VIRTUAL_SITES_RELATIVE
    #define LENNARD_JONES
    #define THOLE
    #define GHOSTS_HAVE_BONDS

.. _Anisotropic non-bonded interactions:

Anisotropic non-bonded interactions
-----------------------------------

..
    .. _Directional Lennard-Jones interaction:

    Directional Lennard-Jones interaction
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. todo::

        Not implemented yet.

    inter LJ-angle

    |image_directional_lj|

    Specifies a 12-10 Lennard-Jones interaction with angular dependence
    between particles of the types ``type1`` and ``type2``. These two particles need two bonded
    partners oriented in a symmetric way. They define an orientation for the
    central particle. The purpose of using bonded partners is to avoid
    dealing with torques, therefore the interaction does *not* need the
    ``ROTATION`` feature. The angular part of the potential minimizes the system
    when the two central beads are oriented along the vector formed by these
    two particles. The shaded beads on the image are virtual particles that
    are formed from the orientation of the bonded partners, connected to the
    central beads. They are used to define angles. The potential is of the
    form

    .. math::

       U(r_{ik},\theta_{jik},\theta_{ikn})=
         \epsilon\left[5\left(\frac{\sigma}r\right)^{12} -
           6\left(\frac{\sigma}{r}\right)^{10}\right]
         \cos^2\theta_{jik}\cos^2\theta_{ikn},

    where :math:`r_{ik}` is the distance between the two central beads, and
    each angle defines the orientation between the direction of a central
    bead (determined from the two bonded partners) and the vector
    :math:`\mathbf{r_{ik}}`. Note that the potential is turned off if one of
    the angle is more than :math:`\pi/2`. This way we don't end up creating
    a minimum for an anti-parallel configuration.

    Unfortunately, the bonded partners are not sought dynamically. One has
    to keep track of the relative positions of the particle IDs. This can be
    done by setting the parameters , , , and . Say the first bead has
    particle ID , then one should set the simulation such as its two bonded
    partners have particle IDs and , respectively. On a linear chain, for
    example, one would typically have and such that the central bead and its
    two bonded partners have position IDs , , and , respectively. This is
    surely not optimized, but once the simulation is set correctly the
    algorithm is very fast.

    It might turn out to be useful in some
    cases to keep force capping during the whole simulation. This is due to
    the very sharp angular dependence for small distance, compared to
    :math:`\sigma`. Two beads might come very close to each other while
    having unfavorable angles such that the interaction is turned off. Then
    a change in the angle might suddenly turn on the interaction and the
    system will blow up (the potential is so steep that one would need
    extremely small time steps to deal with it, which is not very clever for
    such rare events).

    For instance, when modeling hydrogen bonds (N-H...O=C), one can avoid
    simulating hydrogens and oxygens by using this potential. This comes
    down to implementing a HBond potential between N and C atoms.

    The four other optional
    parameters (, , , ) describe a different interaction strength for a
    subset of the simulation box. The box is divided through the plane in
    two different regions: region 1 which creates an interaction with
    strength , region 2 with interaction strength . The 2nd region is
    defined by its -midplane , its total thickness , and the interface width
    . Therefore, the interaction strength is everywhere except for the
    region of the box :math:`z_0-\delta z/2<z<z_0+\delta z/2`. The interface
    width smoothly interpolates between the two regions to avoid
    discontinuities. As an example, one can think of modeling hydrogen bonds
    in two different environments: water, where the interaction is rather
    weak, and in a lipid bilayer, where it is comparatively stronger.

.. _Gay-Berne interaction:

Gay-Berne interaction
~~~~~~~~~~~~~~~~~~~~~

The interface for a Gay-Berne interaction is provided by the :class:`espressomd.interactions.GayBerneInteraction` class. Interaction parameters can be set via::

    system.non_bonded_inter[type1, type2].gay_berne.set_params(**kwargs)

This defines a Gay-Berne potential for prolate and oblate particles
between particles types ``type1`` and ``type2``. The Gay-Berne potential is an
anisotropic version of the classic Lennard-Jones potential, with
orientational dependence of the range :math:`\sigma_0` and the well-depth :math:`\epsilon_0`.

Assume two particles with orientations given by the unit vectors
:math:`\mathbf{\hat{u}}_i` and :math:`\mathbf{\hat{u}}_j` and
intermolecular vector :math:`\mathbf{r} = r\mathbf{\hat{r}}`. If
:math:`r<r_\mathrm{cut}`, then the interaction between these two
particles is given by

.. math::

   V(\mathbf{r}_{ij}, \mathbf{\hat{u}}_i, \mathbf{\hat{u}}_j) = 4
     \epsilon(\mathbf{\hat{r}}_{ij}, \mathbf{\hat{u}}_i,
     \mathbf{\hat{u}}_j) \left( \tilde{r}_{ij}^{-12}-\tilde{r}_{ij}^{-6}
     \right),

otherwise :math:`V(r)=0`. The reduced radius is

.. math::

   \tilde{r}=\frac{r - \sigma(\mathbf{\hat{r}},
       \mathbf{\hat{u}}_i, \mathbf{\hat{u}}_j)+\sigma_0}{\sigma_0},

where

.. math::

   \sigma( \mathbf{\hat{r}}, \mathbf{\hat{u}}_i,
     \mathbf{\hat{u}}_j) = \sigma_{0} \left\{ 1 - \frac{1}{2} \chi \left[
         \frac{ \left( \mathbf{\hat{r}} \cdot \mathbf{\hat{u}}_i +
             \mathbf{\hat{r}} \cdot \mathbf{\hat{u}}_j \right)^{2} }
         {1 + \chi \mathbf{\hat{u}}_i \cdot \mathbf{\hat{u}}_j } +
         \frac{ \left( \mathbf{\hat{r}} \cdot \mathbf{\hat{u}}_i -
             \mathbf{\hat{r}} \cdot \mathbf{\hat{u}}_j \right)^{2} }
         {1 - \chi \mathbf{\hat{u}}_i \cdot \mathbf{\hat{u}}_j}
       \right] \right\}^{-\frac{1}{2}}

and

.. math::

   \begin{gathered}
     \epsilon(\mathbf{\hat{r}}, \mathbf{\hat{u}}_i,
     \mathbf{\hat{u}}_j) = \\
     \epsilon_0 \left( 1- \chi^{2}(\mathbf{\hat{u}}_i
       \cdot \mathbf{\hat{u}}_j)^{2} \right)^{-\frac {\nu}{2}} \left[1-\frac
       {\chi'}{2} \left( \frac { (\mathbf{\hat{r}} \cdot
           \mathbf{\hat{u}}_i+ \mathbf{\hat{r}} \cdot
           \mathbf{\hat{u}}_j)^{2}} {1+\chi' \, \mathbf{\hat{u}}_i \cdot
           \mathbf{\hat{u}}_j }+ \frac {(\mathbf{\hat{r}} \cdot
           \mathbf{\hat{u}}_i-\mathbf{\hat{r}} \cdot
           \mathbf{\hat{u}}_j)^{2}} {1-\chi' \, \mathbf{\hat{u}}_i \cdot
           \mathbf{\hat{u}}_j } \right) \right]^{\mu}.\end{gathered}

The parameters :math:`\chi = \left(k_1^{2} - 1\right)/\left(k_1^{2} + 1\right)`
and :math:`\chi' = \left(k_2^{1/\mu} -  1\right)/\left(k_2^{1/\mu} + 1\right)`
are responsible for the degree of anisotropy of the molecular properties. :math:`k_1` is
the molecular elongation, and :math:`k_2` is the ratio of the potential well depths for the
side-by-side and end-to-end configurations. The exponents and are adjustable
parameters of the potential. Several Gay-Berne parametrizations exist, the
original one being :math:`k_1 = 3`, :math:`k_2 = 5`,
:math:`\mu = 2` and :math:`\nu = 1`.

.. |image_directional_lj| image:: figures/hbond.pdf
