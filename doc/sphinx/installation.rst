.. _Getting, compiling and running:

Getting, compiling and running |es| 
===================================

This chapter will describe how to get, compile and run the software.

|es| releases are available as source code packages from the homepage [1]_.
This is where new users should get the code. The code within release
packages is tested and known to run on a number of platforms.
Alternatively, people that want to use the newest features of or that
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

Example
^^^^^^^

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

Example
^^^^^^^

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
---------------------

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

.. _myconifg.hpp\: Activating and deactivating features:

``myconfig.hpp``: Activating and deactivating features
------------------------------------------------------

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

.. _running a simulation script:


Installing python dependencies
------------------------------

There are a few python packages needed to e.g. build the documentation.
To install the required packages as a non-root user execute the following
command in |es| 's source directory:

.. code-block:: bash

    pip install -r requirements.txt --user --upgrade


Running a simulation script
---------------------------

After |es| is successfully build, a simulation script can fired up
by calling the ``pypresso`` python interpreter located in the build
directory::

    $ ./pypresso <SCRIPT>

The ``pypresso`` script is just a wrapper in order to expose our
self built python modules to the systems python interpreter by
modifying the  ``PYTHONPATH``.
Please see the following chapters describing how to actually write
a simulation script for |es|.


Listing the features compiled into ESPResSo
-------------------------------------------
To see, what features were activated in myconfig.hpp, run:::
    ./pypresso
and then in the Python interpreter:
    import espressomd
    print(espressomd.features())



.. [1]
   http://espressomd.org

.. [2]
   http://git.org

.. [3]
   https://github.com/espressomd/espresso

.. [4]
   https://cmake.org/
