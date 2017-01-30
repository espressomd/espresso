.. include:: defs.rst

First steps
===========

Quick installation
------------------

If you have installed the requirements (see section ) in standard
locations, to compile , it is usually enough to execute the following
sequence of two steps in the directory where you have unpacked the
sources:

::

    $ ./configure make

This will compile in a freshly created object path named according to
your CPU and operating system. As you have not yet specified a
configuration, a standard version will be built with the most often used
features. Usually you will want to build another version of with options
better suited for your purpose.

In some cases, when needs to be compiled for several different platforms
or when different versions with different sets of features are required,
it might be useful to execute the commands not in the source directory
itself, but to start ``configure`` from another directory (see section
). Furthermore, many features of can be selectively turned on or off in
the local configuration header (see section ) before starting the
compilation with ``make``.

The shell script ``configure`` prepares the source code for compilation.
It will determine how to use and where to find the different libraries
and tools required by the compilation process, and it will test what
compiler flags are to be used. The script will find out most of these
things automatically. If something is missing, it will complain and give
hints on how to solve the problem. The configuration process can be
controlled with the help of a number of options that are explained in
section .

The command ``make`` will compile the source code. Depending on the
options passed to the program, ``make`` can also be used for a number of
other things:

-  It can install and uninstall the program to some other directories.
   However, normally it is not necessary to actually *install* to run
   it.

-  It can test for correctness.

-  It can build the documentation.

The details of the usage of ``make`` are described in section .

When these steps have successfully completed, can be started with the
command (see section )

Espresso

where is a Tcl script that tells what to do, and has to be written by
the user. You can find some examples in the ``samples`` folder of the
source code directory. If you want to run in parallel, you should have
compiled to use MPI, and need to tell MPI to run in parallel. The actual
invocation is implementation dependent, but in many cases, such as
OpenMPI, you can use

::

    $ mpirun -n Espresso

where is the number of prcessors to be used.

Running 
--------

Python
~~~~~~

is implemented as a Python module. This means that you need to write a
python script for any task you want to perform with . In this chapter,
the basic structure of the interface will be explained. For a practical
introduction, see the tutorials, which are also part of the
distribution. To use , you need to import the espressomd module in your
Python script. To this end, the folder containing the python module
needs to be in the Python search path. The module is located in the
src/python folder under the build directory. A convenient way to run
python with the correct path is to use the pypresso script located in
the build directory.

::

    $ pypresso simulation.py


Python: Basic concepts
----------------------

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
    bjerrum_length=1) system.actors.add(p3m)

The integrator uses by default the velocity verlet algorithm and is
created by the system class. To perform an integration step, execute

::

    system.integrator.run(steps=100)

.. [1]
   http://www.tcl.tk/man/tcl8.5/tutorial/tcltutorial.html
