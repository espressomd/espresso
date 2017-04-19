.. _introduction:

Introduction
============

|es| is a simulation package designed to perform Molecular Dynamics (MD) and
Monte Carlo (MC) simulations. It is meant to be a universal tool for
simulations of a variety of soft matter systems. |es| features a broad range
of interaction potentials which opens up possibilities for performing
simulations using models with different levels of coarse-graining. It
also includes modern and efficient algorithms for treatment of
electrostatics (P3M, MMM-type algorithms, Maggs algorithm, …),
hydrodynamic interactions (DPD, Lattice-Boltzmann), and magnetic
interactions. It is designed to exploit the capabilities of parallel
computational environments. The program is being continuously extended
to keep the pace with current developments both in the algorithms and
software.

The kernel of |es| is written in C with computational efficiency in mind.
Interaction between the user and the simulation engine is provided via a
Tcl scripting interface. This enables setup of arbitrarily complex
systems which users might want to simulate in future, as well as
modifying simulation parameters during runtime.

.. _Guiding principles:

Guiding principles
------------------

|es| is a tool for performing computer simulation and this user guide
describes how to use this tool. However, it should be borne in mind that
being able to operate a tool is not sufficient to obtain physically
meaningful results. It is always the responsibility of the user to
understand the principles behind the model, simulation and analysis
methods he is using. |es| will *not* do that for you!

It is expected that the users of |es| and readers of this user guide have a
thorough understanding of simulation methods and algorithms they are
planning to use. They should have passed a basic course on molecular
simulations or read one of the renown textbooks,
:cite:`frenkel02b`. It is not necessary to understand
everything that is contained in |es|, but it is inevitable to understand all
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
described in section :ref:`generic_lennard_jones_interaction` where the user
can change all available parameters. Where possible, default values are
avoided, providing the user with the possibility of choice.  |es| cannot be
aware whether your particles are representing atoms or billiard balls, so it
cannot check if the chosen parameters make sense and it is the user’s
responsibility to make sure they do.

On the other hand, flexibility of |es| stems from the employment of a
scripting language at the steering level. Apart from the ability to
modify the simulation and system parameters at runtime, many simple
tasks which are not computationally critical can be implemented at this
level, without even touching the C-kernel. For example, simple
problem-specific analysis routines can be implemented in this way and
made interact with the simulation core. Another example of the program’s
flexibility is the possibility to integrate system setup, simulation and
analysis in one single control script. |es| provides commands to create
particles and set up interactions between them. Capping of forces helps
prevent system blow-up when initially some particles are placed on top
of each other. Using the Tcl interface, one can simulate the randomly
set-up system with capped forces, interactively check whether it is safe
to remove the cap and switch on the full interactions and then perform
the actual productive simulation.

.. _Available simulation methods:

Available simulation methods
----------------------------

provides a number of useful methods. The following table shows the
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

If you believe that the status of a certain method is wrong, please
report so to the developers.

.. tabularcolumns:: |c|c|c|

+--------------------------------+------------------------+------------------+
| **Feature**                    | **Development Status** | **Usage Status** |
+================================+========================+==================+
|             **Integrators**, **Thermostats**, **Barostats**                |
+--------------------------------+------------------------+------------------+
| Velocity-Verlet Integrator     | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| Langevin Thermostat            | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| GHMC Thermostat                | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| DPD Thermostat                 | None                   | Good             |
+--------------------------------+------------------------+------------------+
| Isotropic NPT                  | None                   | Single           |
+--------------------------------+------------------------+------------------+
| NEMD                           | None                   | Group            |
+--------------------------------+------------------------+------------------+
| Quarternion Integrator         | None                   | Good             |
+--------------------------------+------------------------+------------------+
|                                **Interactions**                            |
+--------------------------------+------------------------+------------------+
| Short-range Interactions       | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| Directional Lennard-Jones      | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Gay-Berne Interaction          | None                   | Single           |
+--------------------------------+------------------------+------------------+
| Constraints                    | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| Relative Virtual Sites         | Good                   | Good             |
+--------------------------------+------------------------+------------------+
| Center-of-mass Virtual Sites   | None                   | Good             |
+--------------------------------+------------------------+------------------+
| RATTLE Rigid Bonds             | None                   | Group            |
+--------------------------------+------------------------+------------------+
|                              **Coulomb Interaction**                       |
+--------------------------------+------------------------+------------------+
| P3M                            | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| P3M on GPU                     | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Dipolar P3M                    | Group                  | Good             |
+--------------------------------+------------------------+------------------+
| Ewald on GPU                   | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| MMM1D                          | Single                 | Good             |
+--------------------------------+------------------------+------------------+
| MMM2D                          | Single                 | Good             |
+--------------------------------+------------------------+------------------+
| MMM1D on GPU                   | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| ELC                            | Good                   | Good             |
+--------------------------------+------------------------+------------------+
| MEMD                           | Single                 | Group            |
+--------------------------------+------------------------+------------------+
| ICC*                           | Group                  | Group            |
+--------------------------------+------------------------+------------------+
|                         **Hydrodynamic Interaction**                       |
+--------------------------------+------------------------+------------------+
| Lattice-Boltzmann              | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| Lattice-Boltzmann on GPU       | Group                  | Core             |
+--------------------------------+------------------------+------------------+
| DPD                            | None                   | Good             |
+--------------------------------+------------------------+------------------+
| Shan-Chen Multicomponent Fluid | Group                  | Group            |
+--------------------------------+------------------------+------------------+
| Tunable Slip Boundary          | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Stokesian Dynamics             | Single                 | Single           |
+--------------------------------+------------------------+------------------+
|                             **Analysis**                                   |
+--------------------------------+------------------------+------------------+
| uwerr                          | None                   | Good             |
+--------------------------------+------------------------+------------------+
|                              **Input/Output**                              |
+--------------------------------+------------------------+------------------+
| Blockfiles                     | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| VTF output                     | Core                   | Core             |
+--------------------------------+------------------------+------------------+
| VTK output                     | Group                  | Group            |
+--------------------------------+------------------------+------------------+
| PDB output                     | Good                   | Good             |
+--------------------------------+------------------------+------------------+
| Online visulation with VMD     | Good                   | Good             |
+--------------------------------+------------------------+------------------+
|                               **Miscellaneous**                            |
+--------------------------------+------------------------+------------------+
| Grand canonical feature        | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Metadynamics                   | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Parallel Tempering             | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| Electrokinetics                | Group                  | Group            |
+--------------------------------+------------------------+------------------+
| Object-in-fluid                | Group                  | Group            |
+--------------------------------+------------------------+------------------+
| Collision Detection            | Group                  | Group            |
+--------------------------------+------------------------+------------------+
| Catalytic Reactions            | Single                 | Single           |
+--------------------------------+------------------------+------------------+
| mbtools package                | Group                  | Group            |
+--------------------------------+------------------------+------------------+

.. _Basic program structure:

Basic program structure
-----------------------

As already mentioned, |es| consists of two components. The simulation engine
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
(particles, boundary onditions, interactions, ...), control the
simulation, run analysis, and store and load results. The user has at
hand the full readability and functionality of the scripting language.
For instance, it is possible to use the SciPy package for analysis and
PyPlot for plotting. With a certain overhead in efficiency, it can also
be used to reject/accept new configurations in combined MD/MC schemes.
In principle, any parameter which is accessible from the scripting level
can be changed at any moment of runtime. In this way methods like
thermodynamic integration become readily accessible.

The focus of the user guide is documenting the scripting interfacce, its
behaviour and use in the simulation. It only describes certain technical
details of implementation which are necessary for understanding how the
script interface works. Technical documentation of the code and program
structure is contained in the Developers’ guide (see section [sec:dg]).

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
energy scale by Boltzmann’s constant. However, since |es| does not enforce a
particular unit system, we also don’t know the numerical value of the
Boltzmann constant in the current unit system. Therefore, when
specifying the temperature of a thermostat, you actually do not define
the temperature, but the value of the thermal energy :math:`k_B T` in
the current unit system. For example, if you measure energy in units of
:math:`\mathrm{kJ/mol}` and your real temperature should be 300K, then
you need to set the thermostat’s effective temperature to
:math:`k_B 300\, K \mathrm{mol / kJ} = 2.494`.

As long as one remains within the same unit system throughout the whole
|es|-script, there should be no problems.


.. _Requirements:

Requirements
------------

The following libraries and tools are required to be able to compile and
use :

FFTW
    For some algorithms (P:math:`^3`\ M), needs the FFTW library version
    3 or later  [1]_ for Fourier transforms. Again, the header files are
    required.

MPI
    Finally, if you want to use in parallel, you need a working MPI
    environment (that implements the MPI standard version 1.2).


.. Iinstalling Requirements on ubuntu:

Installing Requirements on Ubuntu 16.04 LTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make ESPResSo run on Ubuntu 16.04 LTS, its dependencies can be
installed with:

.. code-block:: bash

    $ sudo apt install build-essential cmake cython python-numpy tcl-dev
    tk-dev libboost-all-dev openmpi-common

Optionally the ccmake utility can be installed for easier configuration:

.. code-block:: bash

    $ sudo apt install cmake-curses-gui


.. _Installing Requirements on Mac OS X:

Installing Requirements on Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make ESPResSo run on Mac OS X 10.9 or higher, its dependencies can be
installed using MacPorts. First, download the installer package
appropriate for your Mac OS X version from
https://www.macports.org/install.php and install it. Then, run the
following commands:

.. code-block:: bash

    $ sudo xcode-select –install sudo xcodebuild -license accept port
    selfupdate port install cmake python27 python27-cython python27-numpy
    tcl tk openmpi-default fftw-3 +openmpi boost +openmpi +python27 port
    select –set cython cython27 port select –set python python27 port select
    –set mpi openmpi-mp-fortran

.. [1]
   http://www.fftw.org/
