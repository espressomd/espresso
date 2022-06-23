.. _Reaction methods:

Reaction methods
================

This chapter describes methods for simulating chemical reaction equilibria
using reactive particles. Chemical species are referred to by an integer value
stored in the particle :attr:`~espressomd.particle_data.ParticleHandle.type`
property. Chemical reactions take place by changing the value in the
:attr:`~espressomd.particle_data.ParticleHandle.type` property via Monte Carlo
moves using the potential energy of the system before and after the reaction
:cite:`turner08a`.

Please keep in mind the following remarks:

* All reaction methods uses Monte Carlo moves which require potential energies.
  Therefore reaction methods require support for energy calculations for all
  active interactions in the simulation. Some algorithms do not support energy
  calculation, e.g. :ref:`OIF <Object-in-fluid>` and
  :ref:`IBM <Immersed Boundary Method for soft elastic objects>`.

* When modeling reactions that do not conserve the number of particles, the
  method has to create or delete particles from the system. This process can
  invalidate particle ids, in which case the particles are no longer numbered
  contiguously. Particle slices returned by ``system.part`` are still iterable,
  but the indices no longer match the particle ids.

* Checkpointing is not supported, since the state of the Mersenne Twister
  random number generator cannot be serialized.

* For improved performance, you can set the type of invalidated particles with
  :meth:`~espressomd.reaction_methods.ReactionAlgorithm.set_non_interacting_type`
  in all reaction method classes.

* Some of the functionality requires particle book-keeping. If your simulation
  script raises runtime errors about "provided particle type X is currently not
  tracked by the system", use :meth:`system.setup_type_map(type_list=[X])
  <espressomd.system.System.setup_type_map>` where ``X`` is the particle
  type to track.

Thermodynamic ensembles
-----------------------

.. _Reaction ensemble:

Reaction ensemble
~~~~~~~~~~~~~~~~~

The reaction ensemble :cite:`smith94c,turner08a` allows to simulate
chemical reactions which can be represented by the general equation:

.. math::

   \mathrm{\nu_1 S_1 +\ \dots\  \nu_l S_l\ \rightleftharpoons\ \nu_m S_m +\ \dots\ \nu_z S_z }
       \label{general-eq}

where :math:`\nu_i` is the stoichiometric coefficient of species
:math:`S_i`. By convention, stoichiometric coefficients of the
species on the left-hand side of the reaction (*reactants*) attain
negative values, and those on the right-hand side (*products*) attain
positive values, so that the reaction can be equivalently written as

.. math::

   \mathrm{\sum_i \nu_i S_i = 0} \,.
       \label{general-eq-sum}


The equilibrium constant of the reaction is then given as

.. math::

   K = \exp(-\Delta_{\mathrm{r}}G^{\ominus} / k_B T)
       \quad\text{with}\quad
       \Delta_{\mathrm{r}}G^{\ominus} = \sum_i \nu_i \mu_i^{\ominus}\,.
       \label{Keq}


Here :math:`k_B` is the Boltzmann constant, :math:`T` is temperature,
:math:`\Delta_{\mathrm{r}}G^{\ominus}` standard Gibbs free energy change
of the reaction, and :math:`\mu_i^{\ominus}` the standard chemical
potential (per particle) of species :math:`i`. Note that thermodynamic equilibrium is
independent of the direction in which we write the reaction. If it is
written with left and right-hand side swapped,
both :math:`\Delta_{\mathrm{r}}G^{\ominus}` and the stoichiometric
coefficients attain opposite signs, and the equilibrium constant attains the inverse value.
Further, note that the equilibrium constant :math:`K` is the
dimensionless *thermodynamic, concentration-based* equilibrium constant,
defined as

.. math::

   K(c^{\ominus}) = (c^{\ominus})^{-\bar\nu} \prod_i (c_i)^{\nu_i}

where :math:`\bar\nu=\sum_i \nu_i`, and :math:`c^{\ominus}` is the reference concentration,
at which the standard chemical potential :math:`\Delta_{\mathrm{r}}G^{\ominus}` was determined.
In practice, this constant is often used with the dimension of :math:`(c^{\ominus})^{\bar\nu}`

.. math::

   K_c(c^{\ominus}) = K(c^{\ominus})\times (c^{\ominus})^{\bar\nu}

A simulation in
the reaction ensemble consists of two types of moves: the *reaction move*
and the *configuration move*. The configuration move changes the configuration
of the system.
In the *forward* reaction, the appropriate number of reactants (given by
:math:`\nu_i`) is removed from the system, and the concomitant number of
products is inserted into the system. In the *backward* reaction,
reactants and products exchange their roles. The acceptance probability
:math:`P^{\xi}` for a move from state :math:`o` to :math:`n` in the reaction
ensemble is given by the criterion :cite:`smith94c`

.. math::

   P^{\xi} = \text{min}\biggl(1,V^{\bar\nu\xi}\Gamma^{\xi}e^{-\beta\Delta E}\prod_{i=1}\frac{N_i^0!}{(N_i^0+\nu_{i}\xi)!}
       \label{eq:Pacc}
       \biggr),

where :math:`\Delta E=E_\mathrm{new}-E_\mathrm{old}` is the change in potential energy,
:math:`V` is the simulation box volume,
:math:`\beta=1/k_\mathrm{B}T` is the Boltzmann factor, and
:math:`\xi` is the extent of reaction, with :math:`\xi=1` for the forward and
:math:`\xi=-1` for the backward direction.

:math:`\Gamma` is proportional to the reaction constant. It is defined as

.. math::

   \Gamma = \prod_i \Bigl(\frac{\left<N_i\right>}{V} \Bigr)^{\bar\nu} = V^{-\bar\nu} \prod_i \left<N_i\right>^{\nu_i} = K_c(c^{\ominus}=1/\sigma^3)

where :math:`\left<N_i\right>/V` is the average number density of particles of type :math:`i`.
Note that the dimension of :math:`\Gamma` is :math:`V^{\bar\nu}`, therefore its
units must be consistent with the units in which |es| measures the box volume,
i.e. :math:`\sigma^3`.

It is often convenient, and in some cases even necessary, that some particles
representing reactants are not removed from or placed at randomly in the system
but their identity is changed to that of the products, or vice versa in the
backward direction.  A typical example is the ionization reaction of weak
polyelectrolytes, where the ionizable groups on the polymer have to remain on
the polymer chain after the reaction.  The replacement rule is that the identity of a given reactant type is
changed to the corresponding product type as long as the corresponding
coefficients allow for it.  Corresponding means having the same position (index) in
the python lists of reactants and products which are used to set up the
reaction.

Multiple reactions can be added to the same instance of the reaction ensemble.

An example script can be found here:

* `Reaction ensemble / constant pH ensemble <https://github.com/espressomd/espresso/blob/python/samples/reaction_methods.py>`_

For a description of the available methods, see :class:`espressomd.reaction_methods.ReactionEnsemble`.

.. _Grand canonical ensemble:

Grand canonical ensemble
~~~~~~~~~~~~~~~~~~~~~~~~

As a special case, all stoichiometric coefficients on one side of the chemical
reaction can be set to zero. Such a reaction creates particles *ex nihilo*, and
is equivalent to exchanging particles with a reservoir. This type of simulation
in the reaction ensemble is equivalent to the grand canonical simulation.
Formally, this can be expressed by the reaction

.. math::

    \mathrm{\emptyset \rightleftharpoons\ \nu_A A  }  \,,

where, if :math:`\nu_A=1`, the reaction constant :math:`\Gamma` defines the chemical potential of species A.
However, if :math:`\nu_A\neq 1`, the statistics of the reaction ensemble becomes
equivalent to the grand canonical only in the limit of large average number of species A in the box.
If the reaction contains more than one product, then the reaction constant
:math:`\Gamma` defines only the sum of their chemical potentials but not the
chemical potential of each product alone.

Since the Reaction Ensemble acceptance transition probability can be
derived from the grand canonical acceptance transition probability, we
can use the reaction ensemble to implement grand canonical simulation
moves. This is done by adding reactions that only have reactants (for the
deletion of particles) or only have products (for the creation of
particles). There exists a one-to-one mapping of the expressions in the
grand canonical transition probabilities and the expressions in the
reaction ensemble transition probabilities.

.. _Constant pH:

Constant pH
~~~~~~~~~~~

As before in the reaction ensemble, one can define multiple reactions (e.g. for an ampholytic system which contains an acid and a base) in one :class:`~espressomd.reaction_methods.ConstantpHEnsemble` instance:

.. code-block:: python

    cpH=reaction_methods.ConstantpHEnsemble(
        temperature=1, exclusion_range=1, seed=77)
    cpH.add_reaction(gamma=K_diss, reactant_types=[0], reactant_coefficients=[1],
                    product_types=[1, 2], product_coefficients=[1, 1],
                    default_charges={0: 0, 1: -1, 2: +1})
    cpH.add_reaction(gamma=1/(10**-14/K_diss), reactant_types=[3], reactant_coefficients=[1], product_types=[0, 2], product_coefficients=[1, 1], default_charges={0:0, 2:1, 3:1} )


An example script can be found here:

* `Reaction ensemble / constant pH ensemble <https://github.com/espressomd/espresso/blob/python/samples/reaction_methods.py>`_

In the constant pH method due to Reed and Reed
:cite:`reed92a` it is possible to set the chemical potential
of :math:`H^{+}` ions, assuming that the simulated system is coupled to an
infinite reservoir. This value is the used to simulate dissociation
equilibrium of acids and bases. Under certain conditions, the constant
pH method can yield equivalent results as the reaction ensemble :cite:`landsgesell17b`. However, it
treats the chemical potential of :math:`H^{+}` ions and their actual
number in the simulation box as independent variables, which can lead to
serious artifacts.
The constant pH method can be used within the reaction ensemble module by
initializing the reactions with the standard commands of the reaction ensemble.

The dissociation constant, which is the input of the constant pH method, is the equilibrium
constant :math:`K_c` for the following reaction:

.. math::

   \mathrm{HA \rightleftharpoons\ H^+ + A^- } \,,

For a description of the available methods, see :class:`espressomd.reaction_methods.ConstantpHEnsemble`.


Widom Insertion (for homogeneous systems)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Widom insertion method measures the change in excess free energy, i.e. the excess chemical potential due to the insertion of a new particle, or a group of particles:

.. math::

   \mu^\mathrm{ex}_B & :=\Delta F^\mathrm{ex} =F^\mathrm{ex}(N_B+1,V,T)-F^\mathrm{ex}(N_B,V,T)\\
   &=-kT \ln \left(\frac{1}{V} \int_V d^3r_{N_B+1} \langle \exp(-\beta \Delta E_\mathrm{pot}) \rangle_{N_B} \right)

For this one has to provide the following reaction to the Widom method:

.. code-block:: python

    type_B=1
    widom = reaction_methods.WidomInsertion(
        temperature=temperature, seed=77)
    widom.add_reaction(reactant_types=[],
    reactant_coefficients=[], product_types=[type_B],
    product_coefficients=[1], default_charges={1: 0})
    widom.calculate_particle_insertion_potential_energy(reaction_id=0)


The call of ``add_reaction`` define the insertion :math:`\mathrm{\emptyset \to type_B}` (which is the 0th defined reaction).
Multiple reactions for the insertions of different types can be added to the same ``WidomInsertion`` instance.
Measuring the excess chemical potential using the insertion method is done by
calling ``widom.calculate_particle_insertion_potential_energy(reaction_id=0)``
multiple times and providing the accumulated sample to
``widom.calculate_excess_chemical_potential(particle_insertion_potential_energy_samples=samples)``.
If another particle insertion is defined, then the excess chemical potential
for this insertion can be measured in a similar fashion by sampling
``widom.calculate_particle_insertion_potential_energy(reaction_id=1)``.
Be aware that the implemented method only works for the canonical ensemble. If the numbers of particles fluctuate (i.e. in a semi grand canonical simulation) one has to adapt the formulas from which the excess chemical potential is calculated! This is not implemented. Also in a isobaric-isothermal simulation (NpT) the corresponding formulas for the excess chemical potentials need to be adapted. This is not implemented.

The implementation can also deal with the simultaneous insertion of multiple particles and can therefore measure the change of excess free energy of multiple particles like e.g.:

.. math::

   \mu^\mathrm{ex, pair}&:=\Delta F^\mathrm{ex, pair}:= F^\mathrm{ex}(N_1+1, N_2+1,V,T)-F^\mathrm{ex}(N_1, N_2 ,V,T)\\
   &=-kT \ln \left(\frac{1}{V^2} \int_V \int_V d^3r_{N_1+1} d^3 r_{N_2+1} \langle \exp(-\beta \Delta E_\mathrm{pot}) \rangle_{N_1, N_2} \right)

Note that the measurement involves three averages: the canonical ensemble average :math:`\langle \cdot \rangle_{N_1, N_2}` and the two averages over the position of particles :math:`N_1+1` and :math:`N_2+1`.
Since the averages over the position of the inserted particles are obtained via brute force sampling of the insertion positions it can be beneficial to have multiple insertion tries on the same configuration of the other particles.

One can measure the change in excess free energy due to the simultaneous insertions of particles of type 1 and 2 and the simultaneous removal of a particle of type 3:

.. math::

   \mu^\mathrm{ex}:=\Delta F^\mathrm{ex, }:= F^\mathrm{ex}(N_1+1, N_2+1, N_3-1,V,T)-F^\mathrm{ex}(N_1, N_2, N_3 ,V,T)

For this one has to provide the following reaction to the Widom method:

.. code-block:: python

    widom.add_reaction(reactant_types=[type_3],
    reactant_coefficients=[1], product_types=[type_1, type_2],
    product_coefficients=[1,1], default_charges={1: 0})
    widom.calculate_particle_insertion_potential_energy(reaction_id=0)

Be aware that in the current implementation, for MC moves which add
and remove particles, the insertion of the new particle always takes
place at the position where the last particle was removed. Be sure
that this is the behavior you want to have. Otherwise implement a new
function ``WidomInsertion::make_reaction_attempt`` in the core.

An example script which demonstrates how to measure the pair excess
chemical potential for inserting an ion pair into a salt solution
can be found here:

* `Widom Insertion <https://github.com/espressomd/espresso/blob/python/samples/widom_insertion.py>`__

For a description of the available methods, see :class:`espressomd.reaction_methods.WidomInsertion`.

Practical considerations
------------------------

.. _Converting tabulated reaction constants to internal units in ESPResSo:

Converting tabulated reaction constants to internal units in |es|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The implementation in |es| requires that the dimension of :math:`\Gamma`
is consistent with the internal unit of volume, :math:`\sigma^3`. The tabulated
values of equilibrium constants for reactions in solution, :math:`K_c`, typically use
:math:`c^{\ominus} = 1\,\mathrm{moldm^{-3}}` as the reference concentration,
and have the dimension of :math:`(c^{\ominus})^{\bar\nu}`. To be used with |es|, the
value of :math:`K_c` has to be converted as

.. math::

   \Gamma = K_c(c^{\ominus} = 1/\sigma^3) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}})
   \Bigl( N_{\mathrm{A}}\bigl(\frac{\sigma}{\mathrm{dm}}\bigr)^3\Bigr)^{\bar\nu}

where :math:`N_{\mathrm{A}}` is the Avogadro number.  For gas-phase reactions,
the pressure-based reaction constant, :math:`K_p` is often used, which can
be converted to :math:`K_c` as

.. math::

   K_p(p^{\ominus}=1\,\mathrm{atm}) = K_c(c^{\ominus} = 1\,\mathrm{moldm^{-3}}) \biggl(\frac{c^{\ominus}RT}{p^{\ominus}}\biggr)^{\bar\nu},

where :math:`p^{\ominus}=1\,\mathrm{atm}` is the standard pressure.
Consider using the python module pint for unit conversion.

.. _Coupling reaction methods to molecular dynamics:

Coupling reaction methods to molecular dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Monte Carlo (MC) sampling of the reaction can  be coupled with a configurational sampling using Molecular Dynamics (MD).
For non-interacting systems this coupling is not an issue, but for interacting systems the insertion of new particles
can lead to instabilities in the MD integration ultimately leading to a crash of the simulation.

This integration instabilities can be avoided by defining a distance around the particles which already exist in the system
where new particles will not be inserted, which is defined by the required keyword ``exclusion_range``.
This prevents big overlaps with the newly inserted particles, avoiding too big forces between particles, which prevents the MD integration from crashing.
The value of the exclusion range does not affect the limiting result and it only affects the convergence and the stability of the integration.  For interacting systems,
it is usually a good practice to choose the exclusion range such that it is comparable to the diameter of the particles.

If particles with significantly different sizes are present, it is desired to define a different exclusion range for each pair of particle types. This can be done by
defining an exclusion radius per particle type by using the optional argument ``exclusion_radius_per_type``. Then, their exclusion range is calculated using
the Lorentz-Berthelot combination rule, *i.e.* ``exclusion_range = exclusion_radius_per_type[particle_type_1] + exclusion_radius_per_type[particle_type_2]``.
If the exclusion radius of one particle type is not defined, the value of the parameter provided in ``exclusion_range`` is used by default.
If the value in ``exclusion_radius_per_type`` is equal to 0, then the exclusion range of that particle type with any other particle is 0.
