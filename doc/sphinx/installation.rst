.. _Getting, compiling and running:

Getting, compiling and running |es| 
===================================

This chapter will describe how to get, compile and run the software.

|es| releases are available as source code packages from the homepage [1]_.
This is where new users should get the code. The code within release packages
is tested and known to run on a number of platforms.
Alternatively, people that want to use the newest features of |es| or that
want to start contributing to the software can instead obtain the
current development code via the version control system software  [2]_
from |es| ’s project page at Github  [3]_. This code might be not as well
tested and documented as the release code; it is recommended to use this
code only if you have already gained some experience in using |es|.

Unlike most other software, no binary distributions of |es| are available,
and the software is usually not installed globally for all users.
Instead, users of |es| should compile the software themselves. The reason for
this is that it is possible to activate and deactivate various features
before compiling the code. Some of these features are not compatible
with each other, and some of the features have a profound impact on the
performance of the code. Therefore it is not possible to build a single
binary that can satisfy all needs. For performance reasons a user
should always activate only those features that are actually needed.
This means, however, that learning how to compile is a necessary evil.
The build system of |es| uses `cmake` [4]_ to compile
software easily on a wide range of platforms.

.. _Requirements:

Requirements
------------

The following tools libraries, including header files, are required to be able
to compile and use ESPResSo:

CMake
    The build system is based on CMake

C++ Compiler
    C++11 capable C++ compiler (e.g., Gcc 4.8.1 or later)

Boost
    A number of advanced C++ features used by ESPResSo is provided by Boost.

FFTW
    For some algorithms (P:math:`^3`\ M), ESPResSo needs the FFTW library
    version 3 or later  [5]_ for Fourier transforms, including header
    files.

MPI
    Because ESPResSo is parallelized with MPI, you need a working MPI
    environment that implements the MPI standard version 1.2.

Python
    ESPResSo's main user interface is via the Python scripting interface. Both, Python 2 and 3 are supported.

Cython
    Cython is used for connecting the C++ core to Python



.. _Installing Requirements on ubuntu:

Installing Requirements on Ubuntu 16.04 LTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make ESPResSo run on Ubuntu 16.04 LTS, its dependencies can be
installed with:

.. code-block:: bash

    sudo apt install build-essential cmake cython python-numpy \
    libboost-all-dev openmpi-common

Optionally the ccmake utility can be installed for easier configuration:

.. code-block:: bash

    $ sudo apt install cmake-curses-gui


.. _Installing Requirements on Mac OS X:

Installing Requirements on Mac OS X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make |es| run on Mac OS X 10.9 or higher, its dependencies can be
installed using MacPorts. First, download the installer package
appropriate for your Mac OS X version from
https://www.macports.org/install.php and install it. Then, run the
following commands:

.. code-block:: bash

    sudo xcode-select --install
    sudo xcodebuild -license accept
    sudo port selfupdate
    sudo port install cmake python27 python27-cython python27-numpy \
    openmpi-default fftw-3 +openmpi boost +openmpi +python27
    sudo port select --set cython cython27
    sudo port select --set python python27
    sudo port select --set mpi openmpi-mp-fortran

Alternatively, you can use Homebrew.

.. code-block:: bash

    sudo xcode-select --install
    sudo xcodebuild -license accept
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install cmake python cython boost boost-mpi fftw
    brew install numpy --without-python3
    ln -s /usr/local/bin/python2 /usr/local/bin/python

Note: If both MacPorts and Homebrew are installed, you will not be able to
run |es|. Therefore, if you have both installed, please uninstall one
or the other by running one of the following two commands:

.. code-block:: bash

    sudo port -f uninstall installed && rm -r /opt/local
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"

Installing python dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a few python packages needed to e.g. build the documentation.
To install the required packages as a non-root user execute the following
command in |es| 's source directory:

.. code-block:: bash

    pip install -r requirements.txt --user --upgrade


.. _quick installation:

Quick installation
------------------

If you have installed the requirements (see section :ref:`Requirements
<requirements>` ) in standard locations, to compile, it is usually enough to
create a build directory and call ``cmake`` and ``make`` (optional steps 
which modify the build process are commented out):

.. code-block:: bash

    mkdir build
    cd build
    #cp myconfig-default.hpp myconfig.hpp # use the default configuration as template
    #nano myconfig.hpp                    # edit to add/remove features as desired
    cmake ..
    #ccmake . // in order to add/remove features like SCAFACOS or CUDA
    make

This will build |es| with a default feature set, namely
:file:`src/core/myconfig-default.hpp`. This file is a ``c++`` header file, 
which defines the features that should be compiled in.
You may want to adjust the feature set to your needs. This can be easily done
by copying the `myconfig-sample.hpp` which has been created in the build 
directory to `myconfig.hpp` and only uncomment the features you want to use in your simulation.

The ``cmake`` command looks for libraries and tools needed by |es|. So |es| 
can only be built if ``cmake`` reports no errors.

The command ``make`` will compile the source code. Depending on the
options passed to the program, ``make`` can also be used for a number of
other things:

*  It can install and uninstall the program to some other directories.
   However, normally it is not necessary to actually *install* to run
   it: ``make install``

*  It can invoke code checks: ``make check`` 

*  It can build this documentation: ``make sphinx``

When these steps have successfully completed, |es| can be started with the
command::

    ./pypresso <SCRIPT>

where is ``<SCRIPT>`` is a ``python`` script which has to
be written by the user. You can find some examples in the :file:`samples`
folder of the source code directory. If you want to run in parallel, you should
have compiled with *Open MPI*, and need to tell MPI to run in parallel. The actual
invocation is implementation dependent, but in many cases, such as
*Open MPI*, you can use

::

    mpirun -n <N> ./pypresso <SCRIPT>

where ``<N>`` is the number of prcessors to be used.


.. _configuring:

Configuring
-----------

.. _myconifg.hpp\: Activating and deactivating features:

``myconfig.hpp``: Activating and deactivating features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|es| has a large number of features that can be compiled into the binary.
However, it is not recommended to actually compile in all possible
features, as this will slow down significantly. Instead, compile in only
the features that are actually required. A strong gain in speed can be
achieved, by disabling all non-bonded interactions except for a single
one, e.g. . For the developers, it is also possible to turn on or off a
number of debugging messages. The features and debug messages can be
controlled via a configuration header file that contains C-preprocessor
declarations. Appendix lists and describes all available features. The
file ``myconfig-sample.hpp`` that configure will generate in the build
directory contains a list of all possible features that can be copied
into your own configuration file. When no configuration header is
provided by the user, a default header, found in
``src/core/myconfig-default.hpp``, will be used that turns on the
default features.

When you distinguish between the build and the source directory, the
configuration header can be put in either of these. Note, however, that
when a configuration header is found in both directories, the one in the
build directory will be used.

By default, the configuration header is called ``myconfig.hpp``.
The configuration header can be used to compile different binary
versions of with a different set of features from the same source
directory. Suppose that you have a source directory ``$srcdir`` and two
build directories ``$builddir1`` and ``$builddir2`` that contain
different configuration headers:

*  ``$builddir1/myconfig.hpp``:

.. code-block:: c

    #define ELECTROSTATICS
    #define LENNARD-JONES

*  ``$builddir2/myconfig.hpp``:

.. code-block:: c

   #define LJCOS

Then you can simply compile two different versions of via::

    cd builddir1
    cmake ..
    make

    cd builddir2
    cmake ..
    make

To see, what features were activated in myconfig.hpp, run:::

    ./pypresso

and then in the Python interpreter::

    import espressomd
    print(espressomd.features())

.. _cmake:



cmake
-----

In order to build the first step is to create a build directory in which
cmake can be executed. In cmake, the *source directory* (that contains
all the source files) is completely separated from the *build directory*
(where the files created by the build process are put). `cmake` is
designed to *not* be executed in the source directory. `cmake` will
determine how to use and where to find the compiler, as well as the
different libraries and tools required by the compilation process. By
having multiple build directories you can build several variants of |es|,
each variant having different activated features, and for as many
platforms as you want.

**Example:**

When the source directory is ``srcdir`` (the files where unpacked to this
directory), then the user can create a build directory ``build`` below that
path by calling ``mkdir srcdir/build``. In the build direcotry `cmake` is to be
executed, followed by a call of make. None of the files in the source directory
is ever modified when by the build process.

.. code-block:: bash

    $ cd build 
    $ cmake .. 
    $ make

Afterwards Espresso can be run via calling ``./pypresso`` from the command
line.

.. _ccmake:

ccmake
------

Optionally and for easier use the curses interface to cmake can be used
to configure |es| interactively.

**Example:**

Alternatively to the previous example instead of , the executable is
called in the build direcotry to configure ESPResSo previous to its
compilation followed by a call of make:

.. code-block:: bash

    $ cd build 
    $ ccmake .. 
    $ make

Fig. :ref:`ccmake-figure` shows the interactive ccmake UI.

.. _ccmake-figure:

.. figure:: figures/ccmake-example.png
   :alt: ccmake interface
   :width: 70.0%
   :align: center

   ccmake interface


.. _Options and Variables:

Options and Variables
~~~~~~~~~~~~~~~~~~~~~

The behaviour of |es| can be controlled by the means of options and variables
in the CMakeLists.txt file. Also options are defined there. The following
options are available:

* WITH\_CUDA: Build with GPU support

* WITH\_HDF5: Build with HDF5

* WITH\_TESTS: Enable tests

* WITH\_SCAFACOS: Build with Scafacos support

* WITH\_VALGRIND\_INSTRUMENTATION: Build with valgrind instrumentation
  markers

When the value in the CMakeLists.txt file is set to ON the corresponding
option is created if the value of the opition is set to OFF the
corresponding option is not created. These options can also be modified
by calling cmake with the command line argument ``-D``::

    cmake -D WITH_HDF5=OFF srcdir

In the rare event when working with cmake and you want to have a totally
clean build (for example because you switched the compiler), remove the
build directory and create a new one.



.. _make\: Compiling, testing and installing:

``make``: Compiling, testing and installing 
--------------------------------------------

The command ``make`` is mainly used to compile the source code, but it
can do a number of other things. The generic syntax of the ``make``
command is:

.. code-block:: bash

    $ make [options] [target] [variable=value]

When no target is given, the target ``all`` is used. The following
targets are available:

``all``
    Compiles the complete source code. The variable can be used to
    specify the name of the configuration header to be used.

``check``
    Runs the testsuite. By default, all available tests will be run on
    1, 2, 3, 4, 6, or 8 processors.
    
``clean``
    Deletes all files that were created during the compilation.

``install``
    Install |es|. 
    Use ``make DESTDIR=/home/john install`` to install to a 
    specific directory.

``doxygen``
    Creates the Doxygen code documentation in the ``doc/doxygen``
    subdirectory.

``sphinx``
    Creates the `sphinx` code documentation in the ``doc/sphinx``
    subdirectory.

``tutorials``
    Creates the tutorials in the ``doc/tutorials`` subdirectory.

``doc``
    Creates all documentation in the ``doc`` subdirectory (only when
    using the development sources).

A number of options are available when calling ``make``. The most
interesting option is probably ``-j num_jobs``, which can be used for
parallel compilation on computers that have more than one CPU or core.
*num\_jobs* specifies the maximal number of jobs that will be run.
Setting *num\_jobs* to the number of available processors speeds up the
compilation process significantly.


.. _Running es:

Running |es|
------------

|es| is implemented as a Python module. This means that you need to write a
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

    ./pypresso simulation.py

The ``pypresso`` script is just a wrapper in order to expose our
self built python modules to the systems python interpreter by
modifying the  ``PYTHONPATH``.
Please see the following chapters describing how to actually write
a simulation script for |es|.


Basic python simulation script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _Debugging es:

Debugging |es|
--------------

Exceptional situations occur in every program.  If |es| crashes with a
segmentation fault that means that there was a memory fault in the
simulation core which requires running the program in a debugger.  The
`pypresso` executable file is acutally not a program but a script
which sets the Python path appropriately and starts the Python
interpreter with your arguments.  Thus it is not possible to directly
run `pypresso` in a debugger.  However, we provide some useful
commandline options for the most common tools.

.. code-block:: bash

     ./pypresso --tool <args>

where ``--tool`` can be any from the following table.  You can only
use one tool at a time.
  
+---------------------+----------------------------------------------+
| Tool                | Effect                                       |
+=====================+==============================================+
| ``--gdb``           | ``gdb --args python <args>``                 |
+---------------------+----------------------------------------------+
| ``--lldb``          | ``lldb -- python <args>``                    |
+---------------------+----------------------------------------------+
| ``--valgrind``      | ``valgrind --leak-check=full python <args>`` |
+---------------------+----------------------------------------------+
| ``--cuda-gdb``      | ``cuda-gdb --args python <args>``            |
+---------------------+----------------------------------------------+
| ``--cuda-memcheck`` | ``cuda-memcheck python <args>``              |
+---------------------+----------------------------------------------+


.. [1]
   http://espressomd.org

.. [2]
   http://git.org

.. [3]
   https://github.com/espressomd/espresso

.. [4]
   https://cmake.org/

.. [5]
   http://www.fftw.org/
