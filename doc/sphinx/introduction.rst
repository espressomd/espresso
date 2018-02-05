.. _Introduction:

Introduction
============

ESPResSo is a simulation package designed to perform Molecular Dynamics (MD) and
Monte Carlo (MC) simulations. It is meant to be a universal tool for
simulations of a variety of soft matter systems. It features a broad range
of interaction potentials which opens up possibilities for performing
simulations using models with different levels of coarse-graining. It
also includes modern and efficient algorithms for treatment of
electrostatics (P3M, MMM-type algorithms, Maggs algorithm, …),
hydrodynamic interactions (DPD, Lattice-Boltzmann), and magnetic
interactions. It is designed to exploit the capabilities of parallel
computational environments. The program is being continuously extended
to keep the pace with current developments both in the algorithms and
software.

The kernel of ESPResSo is written in C++ with computational efficiency in mind.
Interaction between the user and the simulation engine is provided via a
Python scripting interface. This enables setup of arbitrarily complex
systems which users might want to simulate in future, as well as
modifying simulation parameters during runtime.

.. _Guiding principles:

Guiding principles
------------------

ESPResSo is a tool for performing computer simulation and this user guide
describes how to use this tool. However, it should be borne in mind that
being able to operate a tool is not sufficient to obtain physically
meaningful results. It is always the responsibility of the user to
understand the principles behind the model, simulation and analysis
methods he or she is using. 

It is expected that the users of ESPResSo and readers of this user guide have a
thorough understanding of simulation methods and algorithms they are
planning to use. They should have passed a basic course on molecular
simulations or read one of the renown textbooks,
:cite:`frenkel02b`. It is not necessary to understand
everything that is contained in ESPResSo, but it is inevitable to understand all
methods that you want to use. Using the program as a black box without
proper understanding of the background will most probably result in
wasted user and computer time with no useful output.

To enable future extensions, the functionality of the program is kept as
general as possible. It is modularized, so that extensions to some parts
of the program (implementing a new potential) can be done by modifying
or adding only few files, leaving most of the code untouched.

To facilitate the understanding and the extensibility, much emphasis is
put on readability of the code. Hard-coded assembler loops are generally
avoided in hope that the overhead in computer time will be more than
compensated for by saving much of the user time while trying to
understand what the code is supposed to do.

Hand-in-hand with the extensibility and readability of the code comes the
flexibility of the whole program. On the one hand, it is provided by the
generalized functionality of its parts, avoiding highly specialized functions.
An example can be the implementation of the Generic Lennard-Jones potential
described in section :ref:`Generic Lennard-Jones interaction` where the user
can change all available parameters. Where possible, default values are
avoided, providing the user with the possibility of choice.  ESPResSo cannot be
aware whether your particles are representing atoms or billiard balls, so it
cannot check if the chosen parameters make sense and it is the user’s
responsibility to make sure they do.

On the other hand, flexibility of ESPResSo stems from the employment of a
scripting language at the steering level. Apart from the ability to
modify the simulation and system parameters at runtime, many simple
tasks which are not computationally critical can be implemented at this
level, without even touching the C++-kernel. For example, simple
problem-specific analysis routines can be implemented in this way and
made interact with the simulation core. Another example of the program’s
flexibility is the possibility to integrate system setup, simulation and
analysis in one single control script. ESPResSo provides commands to create
particles and set up interactions between them. Capping of forces helps
prevent system blow-up when initially some particles are placed on top
of each other. Using the Python interface, one can simulate the randomly
set-up system with capped forces, interactively check whether it is safe
to remove the cap and switch on the full interactions and then perform
the actual productive simulation.

.. _Basic program structure:

Basic program structure
-----------------------

As already mentioned, ESPResSo consists of two components. The simulation engine
is written in C and C++ for the sake of computational efficiency. The
steering or control level is interfaced to the kernel via an interpreter
of Python scripting languages.

The kernel performs all computationally demanding tasks. Before all,
integration of Newton’s equations of motion, including calculation of
energies and forces. It also takes care of internal organization of
data, storing the data about particles, communication between different
processors or cells of the cell-system. The kernel is modularized so
that basic functions are accessed via a set of well-defined lean
interfaces, hiding the details of the complex numerical algorithms.

The scripting interface (Python) is used to setup the system
(particles, boundary conditions, interactions, ...), control the
simulation, run analysis, and store and load results. The user has at
hand the full readability and functionality of the scripting language.
For instance, it is possible to use the SciPy package for analysis and
PyPlot for plotting. With a certain overhead in efficiency, it can also
be used to reject/accept new configurations in combined MD/MC schemes.
In principle, any parameter which is accessible from the scripting level
can be changed at any moment of runtime. In this way methods like
thermodynamic integration become readily accessible.

The focus of the user guide is documenting the scripting interface, its
behavior and use in the simulation. It only describes certain technical
details of implementation which are necessary for understanding how the
script interface works. Technical documentation of the code and program
structure is contained in the Developers’ guide (see section [sec:dg]).

.. _Basic python simulation script:

Basic python simulation script
------------------------------

In this section, a brief overview is given over the most important
components of the Python interface and their usage is illustrated by
short examples. The interface is contained in the espressomd Python
module, which needs to be imported, before anything related can be done.

::

    import espressomd

Access to the simulation system is provided via the System class. As a
first step, an instance of the class needs to be created

::

    system=espressomd.System()

Note that only one instance of the System class can be created, due to
limitations in the simulation core. Properties of the System class are
used to access the parameters concerning the simulation system as a
whole, , the box geometry and the time step

::

    system.box_l =(10.0,10.0,15.0) print system.time_step

The particles in the simulation are accessed via the ParticleList class.
It is used to retrieve individual particles of the simulation as well as
for adding particles. An instance of the class is provided as the part
attribute of the System class. Individual particles can be retrieved by
their numerical id by using angular brackets

::

    p=system.part[0]

It is also possible to loop over all particles

::

    for p in system.part: ...

Particles are added via the add method

::

    p=system.part.add(id=1,pos=(3.0,0.5,1.0),q=1)

An individual particle is represented by an instance of ParticleHandle.
The properties of the particle are implemented as Python properties:

::

    p=system.part[0] p.pos=(0,0,0) print p.id,p.pos system.part[0].q=-1

Properties of several particles can be accessed by using Python ranges

::

    v=system.part[:].v

Interactions between particles fall in three categories:

-  Non-bonded interactions are short-ranged interactions between *all*
   pairs of particles of specified types. An example is the
   Lennard-Jones interaction mimicking overlap repulsion and van der
   Wals attraction.

-  Bonded interactions act only between two specific particles. An
   example is the harmonic bond between adjacent particles in a polymer
   chain.

-  Long-range interactions act between all particles with specific
   properties in the entire system. An example is the coulomb
   interaction.

Non-bonded interactions are represented as subclasses of
:class:`espressomd.interactions.NonBondedInteraction`, e.g.
:class:`espressomd.interactions.LennardJonesInteraction`.
Instances of these classes for a given pair of particle types are accessed via
the non_bonded_inter attribute of the System class. Parameters are set as
follows

::

    system.non_bonded_inter[0,0].lennard_jones.set_params(epsilon=1,sigma=1,cutoff=1.5,shift=“auto”)

Bonded interactions are represented by subclasses of BondedInteraction.
To set up a bonded interaction, first an instance of the appropriate
class is created with the desired parameters. Then, the bonded
interaction is registered with the simulation core. Finally, the bond
can be added to particles using the add_bond()-method of ParticleHandle
with the instance of the bond class and the id of the bond partner
particle.

::

    from espressomd.interactions import HarmonicBond
    harmonic=HarmonicBond(k=1,r_0=1) system.bonded_inter.add(harmonic)
    system.part[0].add_bond((harmonic,1))
    system.part[1].add_bond((harmonic,2))

Long-range interactions are subclasses of Actor. They are used by first
creating an instance of the desired actor and then adding it to the
system. To activate the P3M electrostatics solver, execute

::

    from espressomd.electrostatics import P3M p3m=P3M(accuracy=1E-3,
    prefactor=1) system.actors.add(p3m)

The integrator uses by default the velocity verlet algorithm and is
created by the system class. To perform an integration step, execute

::

    system.integrator.run(steps=100)

.. _Tutorials:

Tutorials
---------

There is a number of tutorials that guide you in creating simulations using a number of different features of |es|. You can find them in the directory ``/doc/tutorials``.

.. _Sample scripts:

Sample scripts
--------------

Several scripts that can serve as usage examples can be found in the directory ``/samples``

* ``billard.py`` 
    A simple billard game, needs the Python ``pypopengl`` module

* ``bonds-tst.py``
   Test script that manually creates and deletes different bonds between particles (see :ref:`Bonded interactions`). This script performs:
  
   * print defined bonded interactions 
   * print bonds on a particle
   * delete bonds by index or name
   * save/load a bond to/from a variable
 

* ``cellsystem_test.py``
    Test script that changes the skin depth parameter.  This should not be seen as a benchmark, but rather as a rough estimate of the effect of the cellsystem.     
    .. todo:: implement the older [tune_cells] call
    .. todo:: add save/load optimal cell parameters from tune_skin()
    

* ``coulomb_debye_hueckel.py``,  ``debye_hueckel.py``
    Charged beads with a WCA interaction are simulated using the screened Debye-Hückel potential. See :ref:`Debye-Hückel potential`


* ``ekboundaries.py``

* ``electrophoresis.py``

* ``h5md.py``

* ``lbf.py``

* ``lj-demo.py``
    Lennard-Jones liquid used for demonstration purposes to showcase |es|.
    Sliders from a MIDI controller can change system variables such as
    temperature and volume. Some thermodynamic observables are analyzed and
    plotted live.

* ``lj_liquid_distribution.py``
    Uses ``analysis.distribution`` (See :ref:`Particle distribution`) to analyze a simple Lennard-Jones liquid.

* ``lj_liquid.py``
    Simple Lennard-Jones particle liquid. Shows the basic features of how to:

    * set up system parameters, particles and interactions.
    * warm up and integrate. 
    * write parameters, configurations and observables to files

* ``lj_liquid_structurefactor.py``
    Uses ``analysis.structure_factor`` (See :ref:`Structure factor`) to analyze a simple Lennard-Jones liquid.


* ``load_bonds.py``,  ``store_bonds.py``
    Uses the Python ``pickle`` module to store and load bond information.

* ``load_checkpoint.py``,  ``save_checkpoint.py``
   Basing usage of the checkpointing feature. Shows how to write/load the state of:   
   * custom user variables
   * non bonded interactions
   * particles
   * P3M paremeters
   * thermostat

* ``load_properties.py``,  ``store_properties.py``
    Uses the Python ``pickle`` module to store and load system information.

* ``MDAnalysisIntegration.py``.
    Shows how to expose configuration to ``MDAnalysis`` at run time. The functions of ``MDAnalysis`` can be used to perform some analysis or 
    convert the frame to other formats (CHARMM, GROMACS, ...)

* ``minimal-charged-particles.py``
   Simple Lennard-Jones particle liquid where the particles are assigned charges. The P3M method is used to calculate electrostatic interactions. 

* ``minimal-diamond.py``

* ``minimal-polymer.py``
   Sets up a single dilute bead-spring polymer. Shows the basic usage of ``create_polymer``.

* ``minimal_random_number_generator.py``

* ``observables_correlators.py``

* ``p3m.py``
   Simple Lennard-Jones particle liquid where the particles are assigned charges. The P3M method is used to calculate electrostatic interactions. 

* ``slice_input.py``
    Uses python array slicing to set and extract various particle properties.

* ``visualization_bonded.py``
    Opengl visualization for bonds.

* ``visualization_openGL.py``
    Sample for an opengl visualization with user-defined keyboard- and timed callbacks.

* ``visualization.py``
    A visualization for mayavi/opengl of the lj-liquid with interactive plotting.

* ``visualization_npt.py``
    Simple test visualization for the NPT ensemble.

* ``visualization_poisseuille.py``
    Visualization for poisseuille flow with Lattice-Boltzmann.

* ``visualization_constraints.py``
    Constraint visualization with opengl with all available constraints (commented out).

* ``visualization_mmm2d.py``
    A visual sample for a constant potential plate capacitor simulated with mmm2d.

.. _On units:

On units
--------

What is probably one of the most confusing subjects for beginners of ESPResSo is,
that ESPResSo does not predefine any units. While most MD programs specify a set
of units, like, for example, that all lengths are measured in Ångström
or nanometers, times are measured in nano- or picoseconds and energies
are measured in :math:`\mathrm{kJ/mol}`, ESPResSo does not do so.

Instead, the length-, time- and energy scales can be freely chosen by
the user. Once these three scales are fixed, all remaining units are
derived from these three basic choices.

The probably most important choice is the length scale. A length of
:math:`1.0` can mean a nanometer, an Ångström, or a kilometer -
depending on the physical system, that the user has in mind when he
writes his ESPResSo-script. When creating particles that are intended to
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

The final choice is the time (or mass) scale. By default, ESPResSo uses a reduced
mass of 1, so that the mass unit is simply the mass of all particles.
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

By activating the feature MASSES, you can specify particle masses in
the chosen unit system.

A special note is due regarding the temperature, which is coupled to the
energy scale by Boltzmann’s constant. However, since ESPResSo does not enforce a
particular unit system, we also don’t know the numerical value of the
Boltzmann constant in the current unit system. Therefore, when
specifying the temperature of a thermostat, you actually do not define
the temperature, but the value of the thermal energy :math:`k_B T` in
the current unit system. For example, if you measure energy in units of
:math:`\mathrm{kJ/mol}` and your real temperature should be 300K, then
you need to set the thermostat’s effective temperature to
:math:`k_B 300\, K \mathrm{mol / kJ} = 2.494`.

As long as one remains within the same unit system throughout the whole
ESPResSo-script, there should be no problems.

.. _Available simulation methods:

Available simulation methods
----------------------------

ESPResSo provides a number of useful methods. The following table shows the
various methods as well as their status. The table distinguishes between
the state of the development of a certain feature and the state of its
use. We distinguish between five levels:

**Core**
    means that the method is part of the core of ESPResSo, and that it is
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


In the "Tested" column, we note whether there is an integration test for the method.

If you believe that the status of a certain method is wrong, please
report so to the developers.

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
| Isotropic NPT                  | None                   | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Quarternion Integrator         | None                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                                **Interactions**                                         |
+--------------------------------+------------------------+------------------+------------+
| Short-range Interactions       | Core                   | Core             | Partial    |
+--------------------------------+------------------------+------------------+------------+
| Constraints                    | Core                   | Core             | Partial    |
+--------------------------------+------------------------+------------------+------------+
| Relative Virtual Sites         | Good                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Center-of-mass Virtual Sites   | None                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| RATTLE Rigid Bonds             | None                   | Group            | No         |
+--------------------------------+------------------------+------------------+------------+
|                              **Coulomb Interaction**                                    |
+--------------------------------+------------------------+------------------+------------+
| P3M                            | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| P3M on GPU                     | Single                 | Single           | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Dipolar P3M                    | Group                  | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| MMM1D                          | Single                 | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| MMM2D                          | Single                 | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| MMM1D on GPU                   | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| ELC                            | Good                   | Good             | No         | 
+--------------------------------+------------------------+------------------+------------+
| ICC*                           | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                         **Hydrodynamic Interaction**                                    |
+--------------------------------+------------------------+------------------+------------+
| Lattice-Boltzmann              | Core                   | Core             | No         |
+--------------------------------+------------------------+------------------+------------+
| Lattice-Boltzmann on GPU       | Group                  | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
|                              **Input/Output**                                           |
+--------------------------------+------------------------+------------------+------------+
| VTF output                     | Core                   | Core             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| VTK output                     | Group                  | Group            | No         |
+--------------------------------+------------------------+------------------+------------+
| PDB output                     | Good                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| Online visualisation (Mayavi)  | Good                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| Online visualisation (OpenGL)  | Good                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
|                               **Miscellaneous**                                         |
+--------------------------------+------------------------+------------------+------------+
| Grand canonical feature        | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Electrokinetics                | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Collision Detection            | Group                  | Group            | Partial    |
+--------------------------------+------------------------+------------------+------------+
| Catalytic Reactions            | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Reaction Ensemble              | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+

+--------------------------------+------------------------+------------------+------------+
| **No Python support**                                                                   |
+--------------------------------+------------------------+------------------+------------+
| GHMC Thermostat                | Single                 | Single           | Yes        |
+--------------------------------+------------------------+------------------+------------+
| DPD Thermostat                 | None                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| NEMD                           | None                   | Group            | No         |
+--------------------------------+------------------------+------------------+------------+
| Directional Lennard-Jones      | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Gay-Berne Interaction          | None                   | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| MEMD                           | Single                 | Group            | Yes        | 
+--------------------------------+------------------------+------------------+------------+
| DPD                            | None                   | Good             | Yes        |
+--------------------------------+------------------------+------------------+------------+
| Shan-Chen Multicomponent Fluid | Group                  | Group            | No         |
+--------------------------------+------------------------+------------------+------------+
| Tunable Slip Boundary          | Single                 | Single           | Yes        |
+--------------------------------+------------------------+------------------+------------+
| uwerr                          | None                   | Good             | yes        | 
+--------------------------------+------------------------+------------------+------------+
| Blockfiles                     | Core                   | Core             | Partial    |
+--------------------------------+------------------------+------------------+------------+
| Online visualisation (VMD)     | Good                   | Good             | No         |
+--------------------------------+------------------------+------------------+------------+
| Metadynamics                   | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Parallel Tempering             | Single                 | Single           | No         |
+--------------------------------+------------------------+------------------+------------+
| Object-in-fluid                | Group                  | Group            | Yes        |
+--------------------------------+------------------------+------------------+------------+

