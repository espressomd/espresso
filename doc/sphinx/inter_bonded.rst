.. _Bonded interactions:

Bonded interactions
===================

Bonded interactions are configured by the 
:class:`espressomd.interactions.BondedInteractions` class, which is
a member of :class:`espressomd.system.System`. Generally, one may use 
the following syntax to activate and assign a bonded interaction::

    system.bonded_inter.add(bond)
    system.part[pid1].add_bond((bond, pid2...))

In general, one instantiates an interaction object *bond* and subsequently passes it 
to :meth:`espressomd.interactions.BondedInteractions.add`. This will enable the
bonded interaction and allows the user to assign bonds between particle ids *pidX*. 
Bonded interactions are identified by either their *bondid* or their appropriate object.

Defining a bond between two particles always involves three steps:
defining the interaction, adding it to the system and applying it to the particles.
To illustrate this, assume that three particles with ids 42, 43 and 12 already exist.
One could for example create FENE bonds (more information about the FENE bond
is provided in subsection :ref:`FENE bond`) between them using::

    fene = FeneBond(k=1, d_r_max=1)
    system.bonded_inter.add(fene)
    system.part[42].add_bond((fene, 43), (fene, 12))
    system.part[12].add_bond((fene, 43))

This will set up a FENE bond between particles 42 and 43, 42 and 12, and 12 and 43.
Note that the *fene* object specifies the type of bond and its parameters,
the specific bonds are stored within the particles. you can find more 
information regarding particle properties in :ref:`Setting up particles`.

.. _Distance dependent bonds:

Distance dependent bonds
------------------------

.. _FENE bond:

FENE bond
~~~~~~~~~

A FENE (finite extension nonlinear elastic) bond can be instantiated via
:class:`espressomd.interactions.FeneBond`::
    
    from espressomd.interactions import FeneBond
    fene = FeneBond(k = <float>, d_r_max = <float>, r_0 = <float>)

This command creates a bond type identifier with a FENE
interaction. The FENE potential

.. math::

   V(r) = -\frac{1}{2} K \Delta r_\mathrm{max}^2\ln \left[ 1 - \left(
         \frac{r-r_0}{\Delta r_\mathrm{max}} \right)^2 \right]

models a rubber-band-like, symmetric interaction between two particles with magnitude 
:math:`K`, maximal stretching length :math:`\Delta r_0` and equilibrium bond length
:math:`r_0`. The bond potential diverges at a particle distance
:math:`r=r_0-\Delta r_\mathrm{max}` and :math:`r=r_0+\Delta r_\mathrm{max}`.

.. _Harmonic bond:

Harmonic bond
~~~~~~~~~~~~~

A harmonic bond can be instantiated via
:class:`espressomd.interactions.HarmonicBond`::
    
    from espressomd.interactions import HarmonicBond
    hb = HarmonicBond(k = <float>, r_0 = <float>, r_cut = <float>)


This creates a bond type identifier with a classical harmonic
potential. It is a symmetric interaction between two particles. With the 
equilibrium length :math:`r_0` and the magnitude :math:`k`. It is given by

.. math:: V(r) = \frac{1}{2} k \left( r - r_0 \right)^2

The third, optional parameter defines a cutoff radius. Whenever a
harmonic bond gets longer than :math:`r_\mathrm{cut}`, the bond will be reported as broken,
and a background error will be raised.

.. _Harmonic Dumbbell Bond:

Harmonic Dumbbell Bond
~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Requires ROTATION feature.


A harmonic bond can be instantiated via
:class:`espressomd.interactions.HarmonicDumbbellBond`::
    
    from espressomd.interactions import HarmonicDumbbellBond
    hdb = HarmonicDumbbellBond(k1 = <float>, k2 = <float>, r_0 = <float>, r_cut = <float>)


This bond is similar to the normal harmonic bond in such a way that it
sets up a harmonic potential, i.e. a spring, between the two particles.
Additionally the orientation of the first particle in the bond will be aligned along
the distance vector between both particles. This alignment can be
controlled by the second harmonic constant :math:`k2`. Keep in mind that orientation will
oscillate around the distance vector and some kind of
friction needs to be present for it to relax.

The roles of the parameters :math:`k1, r_0, r_\mathrm{cut}` are exactly the same as for the
harmonic bond.

..
    .. _Quartic bond:

    Quartic bond
    ~~~~~~~~~~~~

    .. todo::
        Not implemented.


    inter quartic

    This creates a bond type with identificator with a quartic potential.
    The potential is minimal at particle distance :math:`r=R`. It is given
    by

    .. math:: V(r) = \frac{1}{2} K_0 \left( r - R \right)^2 + \frac{1}{4} K_1 \left( r - R \right)^4

    The fourth, optional, parameter defines a cutoff radius. Whenever a
    quartic bond gets longer than , the bond will be reported as broken, and
    a background error will be raised.

.. _Bonded coulomb:

Bonded coulomb
~~~~~~~~~~~~~~

.. note::

    Require ELECTROSTAICS feature.

A pairwise Coulomb interaction can be instantiated via 
:class:`espressomd.interactions.BondedCoulomb`::

    bonded_coulomb = espressomd.interactions.BondedCoulomb(prefactor = 1.0)
    system.bonded_inter.add(bonded_coulomb)
    system.part[0].add_bond((bonded_coulomb, 1))

This creates a bond with a Coulomb pair potential between particles `0` and `1`. 
It is given by

.. math:: V(r) = \frac{\alpha q_1 q_2}{r},

where `q1` and `q2` are the charges of the bound particles and `alpha` is the
Coulomb prefactor. This interaction has no cutoff and acts independently of other
Coulomb interactions.

.. _Subtract P3M short-range bond:

Subtract P3M short-range bond
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    requires the P3M feature.

This bond can be instantiated via
:class:`espressomd.interactions.BondedCoulombP3MSRBond`::
    
    from espressomd.interactions import BondedCoulombP3MSRBond
    subtr_p3m_sr = BondedCoulombP3MSRBond(q1q2 = <float>)

The parameter `q1q2` sets the charge factor of the short-range P3M interaction.
It can differ from the actual partice charges.  This specialized bond can be
used to cancel or add **only the short-range** electrostatic part 
of the P3M solver. A use case is descibed in :ref:`Particle polarizability with
thermalized cold Drude oszillators`.

.. _Subtracted Lennard-Jones bond:

Subtracted Lennard-Jones bond
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`espressomd.interactions.SubtLJ` can be used to exclude certain particle 
pairs from the non-bonded :ref:`Lennard-Jones interaction`:: 

    subtLJ = espressomd.interactions.SubtLJ()
    system.bonded_inter.add(subtLJ)
    system.part[0].add_bond((subt, 1))

This bond subtracts the type-pair specific Lennard-Jones interaction between
the involved particles. This interaction is useful when using other bond
potentials which already include the short-ranged repulsion. This often the
case for force fields or in general tabulated potentials.

.. _Rigid bonds:

Rigid bonds
~~~~~~~~~~~

.. note::

    Requires BOND_CONSTRAINT feature.


A rigid bond can be instantiated via
:class:`espressomd.interactions.RigidBond`::
    
    from espressomd.interactions import RigidBond
    rig = RigidBond(r = <float>, ptol = <float>, vtol = <float> )

To simulate rigid bonds, |es| uses the Rattle Shake algorithm which satisfies
internal constraints for molecular models with internal constraints,
using Lagrange multipliers.:cite:`andersen83a` The constrained bond distance 
is named :math:`r`, the positional tolerance is named :math:`ptol` and the velocity tolerance
is named :math:`vtol`.

.. _Tabulated bond interactions:

Tabulated bond interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    
    required TABULATED feature.


A tabulated bond can be instantiated via
:class:`espressomd.interactions.Tabulated`::
    
    from espressomd.interactions import Tabulated
    tab = Tabulated(type = <str>, min = <min>, max = <max>,
                    energy = <energy>, force = <force>)

This creates a bond type identifier with a two-body bond length, 
three-body angle or four-body dihedral 
tabulated potential. For details of the interpolation, see :ref:`Tabulated interaction`.

The bonded interaction can be based on a distance, a bond angle or a
dihedral angle. This is determined by the ``type`` argument, which can
be one of the strings ``distance``, ``angle`` or ``dihedral``.

.. _Calculation of the force and energy:

Calculation of the force and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The potential is calculated as follows:

-  ``type="distance"``: is a two body interaction
   depending on the distance of two particles. The force acts in the
   direction of the connecting vector between the particles. The bond
   breaks above the tabulated range, but for distances smaller than the
   tabulated range, a linear extrapolation based on the first two
   tabulated force values is used.

-  ``type="angle"``: is a three-body angle
   interaction similar to the bond angle potential.
   It is assumed that the potential is tabulated
   for all angles between 0 and :math:`\pi`, where 0 corresponds to a
   stretched polymer, and just as for the tabulated pair potential, the
   forces are scaled with the inverse length of the connecting vectors.
   The force on the extremities acts perpendicular 
   to the connecting vector
   between the corresponding particle and the center particle, in the plane
   defined by the three particles. The force on the center particle
   :math:`p_2` balances the other two forces.

-  ``type="dihedral"``: tabulates a torsional
   dihedral angle potential. It is assumed
   that the potential is tabulated for all angles between 0 and
   :math:`2\pi`. *This potential is not tested yet! Use on own risk, and
   please report your findings and eventually necessary fixes.*

.. _Virtual bonds:

Virtual bonds
~~~~~~~~~~~~~

A virtual bond can be instantiated via
:class:`espressomd.interactions.Virtual`::
    
    from espressomd.interactions import Virtual
    tab = Virtual()


This creates a virtual bond type identifier for a pair bond
without associated potential or force. It can be used to specify topologies
and for some analysis that rely on bonds, or for bonds that should be
displayed in the visualization.

.. _Bond-angle interactions:

Bond-angle interactions
-----------------------
.. note::
    `Feature BOND_ANGLE required.`

Bond-angle interactions involve three particles forming the angle :math:`\phi`, as shown in the schematic below.

.. _inter_angle:
.. figure:: figures/inter_angle.png
   :alt: Bond-angle interactions
   :align: center
   :height: 12.00cm

This allows for a bond type having an angle dependent potential.
This potential is defined between three particles.
The particle for which the bond is created, is the central particle, and the
angle :math:`\phi` between the vectors from this particle to the two
others determines the interaction.

Similar to other bonded interactions, these are defined for every particle triad and and must be added to a particle (see :attr:`espressomd.particle_data.ParticleHandle.bonds`).
For example, for the schematic with particles ``id=0``, ``1`` and ``2`` the bond was defined using ::

    >>> system.part[1].add_bond((bond_angle, 0, 2))

The parameter ``bond_angle`` is a bond type identifier of three possible bond-angle classes, described below.


:class:`espressomd.interactions.AngleHarmonic`
    A classical harmonic potential of the form: 
    
    .. math:: V(\phi) = \frac{K}{2} \left(\phi - \phi_0\right)^2.

    :math:`K` is the bending constant,
    and the optional parameter :math:`\phi_0` is the equilibrium bond angle in
    radians ranging from 0 to :math:`\pi`.

    If this parameter is not given, it defaults to :math:`\phi_0 = \pi`,
    which corresponds to a stretched conformation.

    Unlike the two other variants, this potential has a kink at
    :math:`\phi=\phi_0+\pi` and accordingly a discontinuity in the
    force, and should therefore be used with caution.

    example ::
        >>> angle_harmonic=AngleHarmonic(bend=1.0, phi0=np.pi)
        >>> system.bonded_inter.add(angle_harmonic)
        >>> system.part[1].add_bond((angle_harmonic, 0, 2))



:class:`espressomd.interactions.AngleCosine`

    Cosine bond angle potential of the form:

    .. math:: V(\phi) = K \left[1 - \cos(\phi - \phi0)\right]

    :math:`K` is the bending constant,
    and the optional parameter :math:`\phi_0` is the equilibrium bond angle in
    radians ranging from 0 to :math:`\pi`.

    If this parameter is not given, it defaults to :math:`\phi_0 = \pi`,
    which corresponds to a stretched conformation.

    Around :math:`\phi_0`, this potential is close to a harmonic one
    (both are :math:`1/2(\phi-\phi_0)^2` in leading order), but it is
    periodic and smooth for all angles :math:`\phi`.

    example ::
        >>> angle_cosine=AngleCosine(bend=1.0, phi0=np.pi)
        >>> system.bonded_inter.add(angle_cosine)
        >>> system.part[1].add_bond((angle_cosine, 0, 2))

:class:`espressomd.interactions.AngleCossquare`

    Cosine square bond angle potential of the form:

    .. math:: V(\phi) = \frac{K}{2} \left[\cos(\phi) - \cos(\phi_0)\right]^2

    This form is used for example in the GROMOS96 force field. The
    potential is :math:`1/8(\phi-\phi_0)^4` around :math:`\phi_0`, and
    therefore much flatter than the two potentials before.

    example ::
        >>> angle_cossquare=AngleCossquare(bend=1.0, phi0=np.pi)
        >>> system.bonded_inter.add(angle_cossquare)
        >>> system.part[1].add_bond((angle_cossquare, 0, 2))


.. _Dihedral interactions:

Dihedral interactions
---------------------

Dihedral interactions are available through the :class:`espressomd.interactions.Dihedral` class.

This creates a bond type with identificator with a dihedral potential, a
four-body-potential. In the following, let the particle for which the
bond is created be particle :math:`p_2`, and the other bond partners
:math:`p_1`, :math:`p_3`, :math:`p_4`, in this order. Then, the
dihedral potential is given by

.. math:: V(\phi) = K\left[1 - \cos(n\phi - p)\right],

where :math:`n` is the multiplicity of the potential (number of minima) and can
take any integer value (typically from 1 to 6), :math:`p` is a phase
parameter and :math:`K` is the bending constant of the potential. :math:`\phi` is
the dihedral angle between the particles defined by the particle
quadrupel :math:`p_1`, :math:`p_2`, :math:`p_3` and :math:`p_4`, the
angle between the planes defined by the particle triples :math:`p_1`,
:math:`p_2` and :math:`p_3` and :math:`p_2`, :math:`p_3` and
:math:`p_4`:

|image_dihedral|

Together with appropriate Lennard-Jones interactions, this potential can
mimic a large number of atomic torsion potentials.

.. |image_dihedral| image:: figures/dihedral-angle.pdf


.. _Thermalized distance bond:

Thermalized distance bond
-------------------------

This bond can be used to apply Langevin thermalization on the centre of mass
and the distance of a particle pair.  Each thermostat can have it's own
temperature and friction coefficient. 

The bond is configured with::

    from espressomd.interactions import ThermalizedBond
    thermalized_bond = ThermalizedBond(temp_com = <float>, gamma_com = <float>, temp_distance = <float>, gamma_distance = <float>, r_cut = <float>)
    system.bonded_inter.add(drude_bond)

The parameters are:

    * temp_com : Temerature of the Langevin thermostat for the COM of the particle pair.
    * gamma_com: Friction coefficient of the Langevin thermostat for the COM of the particle pair.
    * temp_distance: Temerature of the Langevin thermostat for the distance vector of the particle pair.
    * gamma_distance: Friction coefficient of the Langevin thermostat for the distance vector of the particle pair.
    * r_cut:  Specifies maximum distance beyond which the bond is considered broken.

The bond is closely related to simulating :ref:`Particle polarizability with
thermalized cold Drude oszillators`. 

