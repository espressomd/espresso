.. _Introduction:

Introduction
============

|es| is a simulation package designed to perform Molecular Dynamics (MD) and
Monte Carlo (MC) simulations. It is meant to be a universal tool for
simulations of a variety of soft matter systems. It :ref:`features<Features>` a
broad range of interaction potentials which opens up possibilities for
performing simulations using models with different levels of coarse-graining.
It also includes modern and efficient algorithms for treatment of
:ref:`electrostatics` (P3M, MMM-type algorithms, constant potential
simulations, dielectric interfaces, …), hydrodynamic interactions
(:ref:`DPD<Dissipative Particle Dynamics (DPD)>`, :ref:`Lattice-Boltzmann`),
and :ref:`magnetic interactions<Magnetostatics>`, only
to name a few.  It is designed to exploit the capabilities of parallel
computational environments.  The program is being continuously extended to keep
the pace with current developments both in the algorithms and software.

The kernel of |es| is written in C++ with computational efficiency in mind.
Interaction between the user and the simulation engine is provided via a
*Python scripting interface*. This enables setup of arbitrarily complex
systems, with simulation parameters that can be modified at runtime.

.. _Guiding principles:

Guiding principles
------------------

|es| is a tool for performing computer simulation and this user guide describes
how to use this tool. However, it should be borne in mind that being able to
operate a tool is not sufficient to obtain physically meaningful results. It is
always the responsibility of the user to understand the principles behind the
model, simulation and analysis methods he or she is using.

It is expected that the users of |es| and readers of this user guide have a
thorough understanding of simulation methods and algorithms they are planning
to use. They should have passed a basic course on molecular simulations or read
one of the renown textbooks, e.g. :cite:`frenkel02b`. It is not necessary to
understand everything that is contained in |es|, but it is inevitable to
understand all methods that you want to use. Using the program as a black box
without proper understanding of the background will most probably result in
wasted user and computer time with no useful output.

To enable future extensions, the functionality of the program is kept as
general as possible. It is modularized, so that extensions to some parts
of the program (e.g. implementing a new potential) can be done by modifying
or adding only a few files, leaving most of the code untouched.

Much emphasis is put on readability of the code. To cite a few examples,
hard-coded C-style for loops are generally avoided in favor of modern C++
range-based for loops or STL accumulators and algorithms, and output
parameters are often avoided by returning a ``std::tuple``. In addition,
vector algebra can be expressed in few lines of code thanks to the
``Utils::Vector`` class that provides overloads for elementary operations,
the dot product, the cross product and operations with matrices.

Hand-in-hand with the extensibility and readability of the code comes the
flexibility of the whole program. On the one hand, it is provided by the
generalized functionality of its parts, avoiding highly specialized functions.
An example can be the implementation of the Generic Lennard-Jones potential
described in section :ref:`Generic Lennard-Jones interaction` where the user
can change all available parameters. Where possible, default values are
avoided, providing the user with the possibility of choice. |es| cannot be
aware whether your particles are representing atoms or billiard balls, so it
cannot check if the chosen parameters make sense and it is the user's
responsibility to make sure they do. In fact, |es| can be used to play
billiard (see sample script :file:`samples/billiard.py`)!

On the other hand, flexibility of |es| stems from the employment of a scripting
language at the steering level. Apart from the ability to modify the simulation
and system parameters at runtime, many simple tasks which are not
computationally critical can be implemented at this level, without even
touching the C++ kernel. For example, simple problem-specific analysis routines
can be implemented in this way and made to interact with the simulation core.
Another example of the program's flexibility is the possibility to integrate
system setup, simulation and analysis in one single control script. |es|
provides commands to create particles and set up interactions between them.
Capping of forces helps prevent system blow-up when initially some particles
are placed on top of each other. Using the Python interface, one can simulate
the randomly set-up system with capped forces, interactively check whether it
is safe to remove the cap and switch on the full interactions and then perform
the actual productive simulation.

.. _Basic program structure:

Basic program structure
-----------------------

As already mentioned, |es| consists of two components. The simulation engine is
written in C++ for the sake of computational efficiency. The steering or
control level is interfaced to the kernel via an interpreter of Python
scripting languages.

The kernel performs all computationally demanding tasks. Before all,
integration of Newton's equations of motion, including calculation of energies
and forces. It also takes care of internal organization of data, storing the
data about particles, communication between different processors or cells of
the cell-system. The kernel is modularized so that basic functions are accessed
via a set of well-defined lean interfaces, hiding the details of the complex
numerical algorithms.

The scripting interface (Python) is used to setup the system (particles,
boundary conditions, interactions, ...), control the simulation, run analysis,
and store and load results. The user has at hand the full readability and
functionality of the scripting language.  For instance, it is possible to use
the SciPy package for analysis and PyPlot for plotting. With a certain overhead
in efficiency, it can also be used to reject/accept new configurations in
combined MD/MC schemes.  In principle, any parameter which is accessible from
the scripting level can be changed at any moment of runtime. In this way
methods like thermodynamic integration become readily accessible.

The focus of the user guide is documenting the scripting interface, its
behavior and use in the simulation. It only describes certain technical details
of implementation which are necessary for understanding how the script
interface works. Technical documentation of the code and program structure is
contained in the `online wiki <https://github.com/espressomd/espresso/wiki>`_.

.. _Basic python simulation script:

Basic python simulation script
------------------------------

In this section, a brief overview is given over the most important components
of the Python interface. Their usage is illustrated by short examples, which
can be put together to a demo script.

.. rubric:: Imports

As usual, the Python script starts by importing the necessary modules.  The
|es| interface is contained in the :mod:`espressomd` Python module, which needs to be
imported, before anything related can be done. ::

    import espressomd

This should be followed by further necessary imports of the example at hand: ::

    import espressomd.interactions
    import espressomd.electrostatics

.. rubric:: espressomd.System

Access to the simulation system is provided via the :class:`~espressomd.system.System` class. As a
first step, an instance of this class needs to be created. ::

    system = espressomd.System(box_l=[10, 10, 10])

Note that only one instance of the System class can be created due to
limitations in the simulation core. :ref:`Properties of the System
class<Setting global variables>` are used to access the parameters
concerning the simulation system such as box geometry, time step or
:ref:`cell-system<Cell systems>`: ::

    print(f"The box dimensions are {system.box_l}")
    system.time_step = 0.01
    system.cell_system.skin = 0.4

.. rubric:: Particles

The particles in the simulation are accessed via ``system.part``, an instance of the
:class:`~espressomd.particle_data.ParticleList` class. Use the ``add`` method to
:ref:`create new particles<Adding particles>`: ::

    part1 = system.part.add(pos=[1.0, 1.0, 1.0], type=0)
    part2 = system.part.add(pos=[1.0, 1.0, 2.0], type=0)

The individual particles are represented by instances of :class:`~espressomd.particle_data.ParticleHandle` which, as
demonstrated in the example above, can be stored as Python variables (``part1`` and ``part2``).
The properties of the particle are implemented as Python properties and can be accessed and/or modified using
the respective :class:`~espressomd.particle_data.ParticleHandle`: ::

    >>> print(part2.pos)
    [1.0, 1.0, 2.0]
    >>> part2.pos = [0.2, 2.0, 0.0]
    >>> print(part2.pos)
    [0.2, 2.0, 0.0]

It is also possible to :ref:`loop<Iterating over particles and pairs of
particles>` over all particles::

    for p in system.part:
        print(f"Particle pos {p.pos}, type {p.type}")

Internally, each particle is automatically assigned a unique numerical id by |es|.
Note that in principle it is possible to explicitly set this particle id (if not in use already) on particle creation.
Using the id, the respective particle can be accessed from the particle list::

    >>> system.part.add(id=3, pos=[2.1, 1.2, 3.3], type=0)
    >>> system.part.by_id(3).pos = [1.0, 1.0, 2.0]
    >>> print(system.part.by_id(3).pos)
    [1.0, 1.0, 2.0]

For larger simulation setups, explicit handling of numerical ids can quickly
become confusing and is thus error-prone. We therefore highly recommend using
:class:`~espressomd.particle_data.ParticleHandle` instead wherever possible.

:ref:`Properties of all particles<Interacting with groups of particles>`
can be accessed via: ::

    positions = system.part.all().pos

.. rubric:: Interactions

In |es|, interactions between particles usually fall in three categories:

-  :ref:`Non-bonded interactions` are short-ranged interactions between *all*
   pairs of particles of specified types. An example is the
   Lennard-Jones interaction mimicking overlap repulsion and van-der-Waals attraction.

-  :ref:`Bonded interactions` act only between two specific particles. An
   example is the harmonic bond between adjacent particles in a polymer
   chain.

-  Long-range interactions act between all particles with specific
   properties in the entire system. An example is the :ref:`Coulomb
   interaction<Electrostatics>`.

.. rubric:: Non-bonded interaction

Non-bonded interactions are represented as subclasses of
:class:`~espressomd.interactions.NonBondedInteraction`, e.g.
:class:`~espressomd.interactions.LennardJonesInteraction`.
Instances of these classes for a given pair of particle types are accessed via
the non_bonded_inter attribute of the System class. This sets up a Lennard-Jones
interaction between all particles of type 0 with the given parameters: ::

    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=1, sigma=1, cutoff=5.0, shift="auto")

.. rubric:: Bonded interaction

Next, we add a pair of particles with a different type to later add
a :ref:`harmonic bond<Harmonic bond>` between them: ::

    part1 = system.part.add(pos=[7.0, 7.0, 7.0], type=1)
    part2 = system.part.add(pos=[7.0, 7.0, 8.0], type=1)

To set up a bonded interaction, first an instance of the appropriate
class is created with the desired parameters: ::

    harmonic = espressomd.interactions.HarmonicBond(k=1.0, r_0=0.5)

Then, the bonded interaction is registered in the simulation core
by adding the instance to :attr:`~espressomd.system.System.bonded_inter`: ::

    system.bonded_inter.add(harmonic)

Finally, the bond can be added to particles using the :meth:`~espressomd.particle_data.ParticleHandle.add_bond()`
method of :class:`~espressomd.particle_data.ParticleHandle` with the instance of the bond class and the
instance of the partner particle: ::

    part1.add_bond((harmonic, part2))

.. rubric:: Charges

Now we demonstrate how to setup a pair of charged particles treated by the P3M
electrostatics solver. We start by adding the particles: ::

    cation = system.part.add(pos=[4.0, 1.0, 1.0], type=2, q=1.0)
    anion = system.part.add(pos=[6.0, 1.0, 1.0], type=2, q=-1.0)

Long-range interactions and other methods that might be mutually exclusive
are treated as so-called *actors*. They are used by first creating an instance
of the desired actor::

    p3m = espressomd.electrostatics.P3M(accuracy=1e-3, prefactor=1.0)

and then adding it to the system: ::

    print("Tuning p3m ...")
    system.actors.add(p3m)

.. rubric:: Integration

So far we just *added* particles and interactions, but did not propagate the
system. This is done by the *integrator*.  It uses by default the velocity
Verlet algorithm and is already created by the system class. To perform an
integration step, just execute::

    system.integrator.run(1)

Usually, the system is propagated for a number of steps in a loop alongside
with some analysis. In this last snippet, the different energy contributions
of the system are printed: ::

    num_configs = 10
    num_steps = 1000

    for i in range(num_configs):
        system.integrator.run(num_steps)
        energy = system.analysis.energy()
        print(f"System time: {system.time}")
        print(f"Energy of the LJ interaction: {energy['non_bonded']}")
        print(f"Energy of the harmonic bond: {energy['bonded']}")
        print(f"Energy of the Coulomb interaction: {energy['coulomb']}")

.. _Tutorials:

Tutorials
---------

There are a number of tutorials that introduce the use of |es| for different
physical systems. You can also find the tutorials and related scripts in the
directory :file:`/doc/tutorials`.
They are executed with the ``ipypresso`` script.

The following tutorials are available:

* :file:`lennard_jones`: Modelling of a single-component and a two-component Lennard-Jones liquid.
* :file:`visualization`: Using the online visualizers of |es|.
* :file:`error_analysis`: Statistical analysis of simulation results.
* :file:`charged_system`: Modelling of ion condensation around a charged rod.
* :file:`ferrofluid`: Modelling a colloidal suspension of magnetic particles.
* :file:`lattice_boltzmann`: Simulations including hydrodynamic interactions using the lattice-Boltzmann method.
* :file:`raspberry_electrophoresis`: Extended objects in a lattice-Boltzmann fluid, raspberry particles.
* :file:`active_matter`: Modelling of self-propelling particles.
* :file:`electrokinetics`: Modelling electrokinetics together with hydrodynamic interactions.
* :file:`constant_pH`: Modelling the titration of a weak acid using the constant pH method
* :file:`widom_insertion`: Measuring the excess chemical potential of a salt solution using the Widom particle insertion method

The executed notebooks with solutions and plots are periodically deployed
online to the `GitHub Pages <https://espressomd.github.io/tutorials.html>`__.

.. _Sample scripts:

Sample scripts
--------------

Several scripts that can serve as usage examples can be found in the
directory :file:`/samples`.
They are executed with the ``pypresso`` script.

The following samples are available:

.. include:: samples.rst


.. _On units:

On units
--------

What is probably one of the most confusing subjects for beginners of |es| is,
that |es| does not predefine any units. While most MD programs specify a set
of units, like, for example, that all lengths are measured in Ångström
or nanometers, times are measured in nano- or picoseconds and energies
are measured in :math:`\mathrm{kJ/mol}`, |es| does not do so.

Instead, the length-, time- and energy scales can be freely chosen by
the user. Once these three scales are fixed, all remaining units are
derived from these three basic choices.

The probably most important choice is the length scale. A length of
:math:`1.0` can mean a nanometer, an Ångström, or a kilometer -
depending on the physical system, that the user has in mind when he
writes his |es|-script. When creating particles that are intended to
represent a specific type of atoms, one will probably use a length scale
of Ångström. This would mean, that the parameter :math:`\sigma` of the
Lennard-Jones interaction between two atoms would be set to twice the
van-der-Waals radius of the atom in Ångström. Alternatively, one could
set :math:`\sigma` to :math:`2.0` and measure all lengths in multiples
of the van-der-Waals radius. When simulation colloidal particles, which
are usually of micrometer size, one will choose their diameter (or
radius) as basic length scale, which is much larger than the Ångström
scale used in atomistic simulations.

The second choice to be made is the energy scale. One can for example
choose to set the Lennard-Jones parameter :math:`\epsilon` to the energy
in :math:`\mathrm{kJ/mol}`. Then all energies will be measured in that
unit. Alternatively, one can choose to set it to :math:`1.0` and measure
everything in multiples of the van-der-Waals binding energy of the
respective particles.

The final choice is the time (or mass) scale. By default, |es| uses a reduced
mass of 1 for all particles, so that the mass unit is simply the mass of one particle.
Combined with the energy and length scale, this is sufficient to derive
the resulting time scale:

.. math::

    [\mathrm{time}] = [\mathrm{length}]\sqrt{\frac{[\mathrm{mass}]}{[\mathrm{energy}]}}

This means, that if you measure lengths in Ångström, energies in
:math:`k_B T` at 300K and masses in 39.95u, then your time scale is
:math:`\mathring{A} \sqrt{39.95u / k_B T} = 0.40\,\mathrm{ps}`.

On the other hand, if you want a particular time scale, then the mass
scale can be derived from the time, energy and length scales as

.. math::

    [\mathrm{mass}] = [\mathrm{energy}]\frac{[\mathrm{time}]^2}{[\mathrm{length}]^2}.

By activating the feature ``MASS``, you can specify particle masses in
the chosen unit system.

A special note is due regarding the temperature, which is coupled to the
energy scale by Boltzmann's constant. However, since |es| does not enforce a
particular unit system, we also don't know the numerical value of the
Boltzmann constant in the current unit system. Therefore, when
specifying the temperature of a thermostat, you actually do not define
the temperature, but the value of the thermal energy :math:`k_B T` in
the current unit system. For example, if you measure energy in units of
:math:`\mathrm{kJ/mol}` and your real temperature should be 300K, then
you need to set the thermostat's effective temperature to
:math:`k_B 300\, K \mathrm{mol / kJ} = 2.494`.

As long as one remains within the same unit system throughout the whole
|es|-script, there should be no problems.

.. _Available simulation methods:

Available simulation methods
----------------------------

|es| provides a number of useful methods. The following table shows the
various methods as well as their status. The table distinguishes between
the state of the development of a certain feature and the state of its
use. We distinguish between five levels:

**Core**
    means that the method is part of the core of |es|, and that it is
    extensively developed and used by many people.

**Good**
    means that the method is developed and used by independent people
    from different groups.

**Group**
    means that the method is developed and used in one group.

**Single**
    means that the method is developed and used by one person only.

**None**
    means that the method is developed and used by nobody.

**Experimental**
    means that the method might have side effects.

In the "Tested" column, we note whether there is an integration test for the method.

If you believe that the status of a certain method is wrong, please
report so to the developers using the instructions in :ref:`Contributing`.

.. tabularcolumns:: |l|c|c|c|

+--------------------------------+------------------------+------------------+------------+
| **Feature**                    | **Development Status** | **Usage Status** | **Tested** |
+================================+========================+==================+============+
|             **Integrators**, **Thermostats**, **Barostats**                             |
+--------------------------------+------------------------+------------------+------------+
| Velocity-Verlet Integrator     | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Langevin Thermostat            | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Isotropic NpT                  | Experimental           | None             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Quaternion Integrator          | Core                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Stokesian Dynamics             | Single                 | None             | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                                **Interactions**                                         |
+--------------------------------+------------------------+------------------+------------+
| Short-range Interactions       | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Constraints                    | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Relative Virtual Sites         | Good                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| RATTLE Rigid Bonds             | Single                 | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Gay--Berne Interaction         | Experimental           | Experimental     | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                              **Coulomb Interaction**                                    |
+--------------------------------+------------------------+------------------+------------+
| P3M                            | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| P3M on GPU                     | Single                 | Single           | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Dipolar P3M                    | Group                  | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| MMM1D                          | Single                 | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| MMM1D on GPU                   | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| ELC                            | Good                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| ICC*                           | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                         **Hydrodynamic Interaction**                                    |
+--------------------------------+------------------------+------------------+------------+
| Lattice-Boltzmann              | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Lattice-Boltzmann on GPU       | Group                  | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                              **Input/Output**                                           |
+--------------------------------+------------------------+------------------+------------+
| VTF output                     | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| VTK output                     | Group                  | Group            | No         |
+--------------------------------+------------------------+------------------+------------+
| Checkpointing                  | Experimental           | Experimental     | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                              **Visualization**                                          |
+--------------------------------+------------------------+------------------+------------+
| OpenGL visualizer              | Good                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
|                               **Miscellaneous**                                         |
+--------------------------------+------------------------+------------------+------------+
| Electrokinetics                | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Collision Detection            | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Reaction Ensemble              | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Constant pH Method             | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Object-in-fluid                | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Immersed boundary method       | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| DPD                            | Single                 | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| DPD Thermostat                 | Single                 | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+



.. _Software releases:

Software releases
-----------------

|es| releases use the following versioning scheme: ``major.minor.patch_level``.
New features are introduced in major and minor releases, while bugfix releases
only patch bugs without adding or removing features. Since the ``patch_level``
doesn't affect the capabilities of the software, it's common to refer to
releases simply as ``major.minor``.

New users should always choose the latest release. When opting for an
older release, we recommend using the latest bugfix release from that
line (for example 4.0.2 instead of 4.0), unless you need to capture the
behavior of bugs for reproducibility reasons. When filing bug reports
or citing |es|, the version should always be mentioned. See
our policies on :ref:`bug reports <Contributing>` and
:ref:`citing the software <How to cite ESPResSo>` for more details.

Releases from 4.0 onward can be found on
`GitHub <https://github.com/espressomd/espresso/releases>`_.
Older releases from 2.1 to 3.3 can be found in
`GNU Savannah <https://download.savannah.gnu.org/releases/espressomd/>`_.
See our policy on :ref:`API backward compatibility
<Intended interface compatibility between ESPResSo versions>`
if you need more details.

.. _Release workflow:

Release workflow
^^^^^^^^^^^^^^^^

Major and minor releases are branched from the development branch ``python``.
When a version ``X.Y.0`` is released, the ``python`` branch is copied
to a new branch named ``X.Y``, at which point the ``python`` branch is ready
to accept contributions for the ``X.Y+1.0`` release. The ``X.Y`` branch
still gets bugfix releases ``X.Y.1``, ``X.Y.2``, ..., for several months.

`GitHub milestones <https://github.com/espressomd/espresso/milestones>`_
track the progress of each release. They can give you an idea of the changes
in future releases, although it's more convenient to follow the live release
notes in the `wiki <https://github.com/espressomd/espresso/wiki>`_ (listed
under "Planned releases" in the side bar). These notes are updated monthly.
Most users will only be interested in the live release notes of the
planned bugfix release for the version of |es| they're using.

If you're actively developing code for |es|, you might also be interested in
the summaries of the `ESPResSo meetings
<https://github.com/espressomd/espresso/wiki/Offline-Espresso-meeting>`_,
where the core team discusses plans for future releases and feature freezes.

.. _Intended interface compatibility between ESPResSo versions:

Intended interface compatibility between |es| versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With regards to the stability of the Python interface, we have the following
guidelines:

* ``patch_level``: The Python interface will not change if only the
  ``patch_level`` number is different. Example: 4.0.0 :math:`\to` 4.0.1.

* ``minor``: There will be no silent interface changes between two versions
  with different minor numbers, i.e. a simulation script will not silently
  produce different results with the new version. The interface can, however,
  be extended. In important cases, the interface can change in such a way
  that using the old interface produces a clear error message and the
  simulation is terminated. Example: 4.0.2 :math:`\to` 4.1.0.

* ``major``: No guarantees are made for a transition between major versions.
  Example: 4.1.2 :math:`\to` 5.0.

* No guarantees are made with regards to the development branch on GitHub.

* No guarantees are made with respect to the C++ bindings in the simulation core.

.. _How to cite ESPResSo:

How to cite |es|
^^^^^^^^^^^^^^^^

Please cite :cite:t:`weik19a` (BibTeX key ``weik19a`` in :file:`doc/bibliography.bib`)
for |es| 4.0 and later, or both :cite:t:`arnold13a` and :cite:t:`limbach06a`
(BibTeX keys ``arnold13a`` and ``limbach06a`` in :file:`doc/bibliography.bib`)
for |es| 2.0 to 3.3. To find the version number, use the following command:

.. code-block:: bash

    ./pypresso -c "import espressomd.version;print(espressomd.version.friendly())"

A number of algorithms in |es| are fairly advanced and unique to |es|.
The authors of these contributions kindly ask you to cite the relevant
publications, using the BibTeX entries indicated in this user guide.

A complete citation would look like this:

    Simulations were carried out with ESPResSo 4.2[24] using the ICC\*
    algorithm[25].

    | ____________

    | [24] F. Weik, R. Weeber, K. Szuttor *et al.* ESPResSo 4.0 -- an
      extensible software package for simulating soft matter systems.
      *Eur. Phys. J. Spec. Top.* **227**, 1789--1816 (2019).
      doi:\ `10.1140/epjst/e2019-800186-9 <https://doi.org/10.1140/epjst/e2019-800186-9>`_.
    | [25] C. Tyagi, M. Süzen, M. Sega *et al.* An iterative, fast,
      linear-scaling method for computing induced charges on arbitrary
      dielectric boundaries. *J. Chem. Phys.* **132**, 154112 (2010).
      doi:\ `10.1063/1.3376011 <https://doi.org/10.1063/1.3376011>`_.

You may also provide the patch level, when relevant. If you developed code
for |es| and made it available in a publicly accessible repository, you
should consider providing the corresponding URL, for example in a footnote:

    The method was implemented for ESPResSo 4.2.1[24] and the source code is
    available online\ :superscript:`note 1`.

    | ____________

    | :superscript:`note 1` https://github.com/username/espresso/tree/implemented-algorithm

    | [24] F. Weik, R. Weeber, K. Szuttor *et al.* ESPResSo 4.0 -- an
      extensible software package for simulating soft matter systems.
      *Eur. Phys. J. Spec. Top.* **227**, 1789--1816 (2019).
      doi:\ `10.1140/epjst/e2019-800186-9 <https://doi.org/10.1140/epjst/e2019-800186-9>`_.

