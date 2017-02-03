.. _getting, compiling and running:

Getting, compiling and running |es| 
===================================

This chapter will describe how to get, compile and run the software.

releases are available as source code packages from the home page [1]_.
This is where new users should get the code. The code within release
packages is tested and known to run on a number of platforms.
Alternatively, people that want to use the newest features of or that
want to start contributing to the software can instead obtain the
current development code via the version control system software  [2]_
from ’s project page at Github  [3]_. This code might be not as well
tested and documented as the release code; it is recommended to use this
code only if you have already gained some experience in using .

Unlike most other software, no binary distributions of are available,
and the software is usually not installed globally for all users.
Instead, users of should compile the software themselves. The reason for
this is that it is possible to activate and deactivate various features
before compiling the code. Some of these features are not compatible
with each other, and some of the features have a profound impact on the
performance of the code. Therefore it is not possible to build a single
binary that can satisfy all needs. F or performance reasons a user
should always activate only those features that are actually needed.
This means, however, that learning how to compile is a necessary evil.
The build system of uses either the GNU autotools or cmake to compile
software easily on a wide range of platforms.

cmake
-----

In order to build the first step is to create a build directory in which
cmake can be executed. In cmake, the *source directory* (that contains
all the source files) is completely separated from the *build directory*
(where the files created by the build process are put). cmake is
designed to not be executed in the source directory. Cmake will
determine how to use and where to find the compiler, as well as the
different libraries and tools required by the compilation process. By
having multiple build directories you can build several variants of ,
each variant having different activated features, and for as many
platforms as you want.

Example
'''''''

When the source directory is (the files where unpacked to this
directory), then the user can create a build directory below that path
by calling the . In the build direcotry cmake is to be executed,
followed by a call of make. None of the files in the source directory is
ever modified when by the build process.

cd build cmake srcdir make Espresso

Afterwards Espresso can be run via calling Espresso from the command
line.

ccmake
~~~~~~

Optionally and for easier use the curses interface to cmake can be used
to configure ESPResSo interactively.

Example
'''''''

Alternatively to the previous example instead of , the executable is
called in the build direcotry to configure ESPResSo previous to its
compilation followed by a call of make:

cd build ccmake srcdir make Espresso

Fig. [fig:ccmake] shows the interactive ccmake UI.

.. figure:: figures/ccmake-example.png
   :alt: ccmake interface
   :width: 70.0%

   ccmake interface

Options and Variables
~~~~~~~~~~~~~~~~~~~~~

The behaviour of can be controlled by the means of options and variables
in the CMakeLists.txt file. Also options are defined there. The
following options are available:

WITH\_PYTHON: Build python interface

WITH\_TCL: Build tcl interface

WITH\_CUDA: Build with GPU support

WITH\_HDF5: Build with HDF5

WITH\_TESTS: Enable tests

WITH\_SCAFACOS: Build with Scafacos support

WITH\_VALGRIND\_INSTRUMENTATION: Build with valgrind instrumentation
markers

When the value in the CMakeLists.txt file is set to ON the corresponding
option is created if the value of the opition is set to OFF the
corresponding option is not created. These options can also be modified
by calling cmake with the command line argument -D:

cmake -D WITH\_TCL=OFF srcdir

In the rare event when working with cmake and you want to have a totally
clean build (for example because you switched the compiler), remove the
build directory and create a new one.

``make``: Compiling, testing and installing 
--------------------------------------------

The command ``make`` is mainly used to compile the source code, but it
can do a number of other things. The generic syntax of the ``make``
command is:

make [] [...] [=]

When no target is given, the target ``all`` is used. The following
targets are available:

``all``
    Compiles the complete source code. The variable can be used to
    specify the name of the configuration header to be used.

``check``
    | Runs the testsuite. By default, all available tests will be run on
      1, 2, 3, 4, 6, or 8 processors. Which tests are run can be
      controlled by means of the variable ``tests``, which processor
      numbers are to be used can be controlled via the variable
      ``processors``. Note that depending on your MPI installation, MPI
      jobs can only be run in the queueing system, so that will not run
      from the command line. In that case, you may not be able to run
      the testsuite, or you have to directly submit the testsuite script
      ``testsuite/test.sh`` to the queueing system.
    | **Example:** ``make check tests="madelung.tcl" processors="1 2"``
    | will run the test ``madlung.tcl`` on one and two processors.

``clean``
    Deletes all files that were created during the compilation.

``mostlyclean``
    Deletes most files that were created during the compilation. Will
    keep for example the built doxygen documentation and the binary.

``dist``
    | Creates a ``.tar.gz``-file of the sources. This will include all
      source files as they currently are in the source directory, it
      will include local changes. This is useful to give your version of
      to other people. The variable ``extra`` can be used to specify
      additional files and directories that are to be included in the
      archive file.
    | **Example:** ``make dist extra="myconfig.hpp internal"``
    | will create the archive file and include the file ``myconfig.hpp``
      and the directory ``internal`` with all files and subdirectories.

``install``
    | Install . The variables ``prefix`` and ``exec-prefix`` can be used
      to specify the installation directories, otherwise the defaults
      defined by the ``configure`` script are used. ``prefix`` sets the
      prefix where all files are to be installed, ``exec-prefix`` sets
      the prefix where the executable files are to be installed and is
      required only when there is an architecture-specific directory.
    | **Example:** ``make install prefix=/usr/local``
    | will install all files below ``/usr/local``.

``ug  ``
    Creates the User guide in the ``doc/ug`` subdirectory (only when
    using the development sources).

``dg  ``
    Creates the Developers’ guide in the ``doc/dg`` subdirectory (only
    when using the development sources).

``doxygen  ``
    Creates the Doxygen code documentation in the ``doc/doxygen``
    subdirectory.

``tutorials  ``
    Creates the tutorials in the ``doc/tutorials`` subdirectory.

``doc ``
    Creates all documentation in the ``doc`` subdirectory (only when
    using the development sources).

A number of options are available when calling ``make``. The most
interesting option is probably ``-j num_jobs``, which can be used for
parallel compilation on computers that have more than one CPU or core.
*num\_jobs* specifies the maximal number of jobs that will be run.
Setting *num\_jobs* to the number of available processors speeds up the
compilation process significantly.

TCL: Running 
-------------

When is found in your path, it can be run via

Espresso [ []]

When is called without any arguments, it is started in the interactive
mode, where new commands can be entered on the command line. When the
name of a *tcl\_script* is given, the script is executed. Any further
arguments are passed to the script.

If you want to run in parallel using MPI, the actual invocation depends
on your MPI implementation. In many cases, OpenMPI, the command will be

mpiexec -n Espresso [ []]

where denotes the number of MPI nodes to be used. However, note that
depending on your MPI installation, MPI jobs can only be run in a
queueing system, so that will not run from the command line. Also, older
installations sometimes require “-np” instead of “-n” or “mpirun”
instead of “mpiexec”.

``myconfig.hpp``: Activating and deactivating features
------------------------------------------------------

has a large number of features that can be compiled into the binary.
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

By default, the configuration header is called ``myconfig.hpp``. The
name of the configuration header can be changed either when the
``configure``-script is called via the variable (see section ), or when
``make`` is called with the setting (see section ).

The configuration header can be used to compile different binary
versions of with a different set of features from the same source
directory. Suppose that you have a source directory ``$srcdir`` and two
build directories ``$builddir1`` and ``$builddir2`` that contain
different configuration headers:

-  ``$builddir1/myconfig.hpp``:

   #define ELECTROSTATICS #define LENNARD-JONES

-  ``$builddir2/myconfig.hpp``:

   #define LJCOS

Then you can simply compile two different versions of via::

    cd builddir1
    srcdir/configure
    make

    cd builddir2
    srcdir/configure
    make

.. [1]
   http://espressomd.org

.. [2]
   http://git.org

.. [3]
   https://github.com/espressomd/espresso
