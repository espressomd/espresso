.. _Installation:

Installation
============

This chapter will describe how to get, compile and run the software.

|es| releases are available as source code packages from the homepage [1]_.
This is where new users should get the code. The code within release packages
is tested and known to run on a number of platforms.
Alternatively, people that want to use the newest features of |es| or that
want to start contributing to the software can instead obtain the
current development code via the version control system software  [2]_
from |es|'s project page at Github  [3]_. This code might be not as well
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
The build system of |es| uses ``cmake`` [4]_ to compile
software easily on a wide range of platforms.

.. _Requirements:

Requirements
------------

The following tools libraries, including header files, are required to be able
to compile and use |es|:

CMake
    The build system is based on CMake

C++ Compiler
    C++14 capable C++ compiler (e.g., gcc 5 or later)

Boost
    A number of advanced C++ features used by |es| are provided by Boost.
    We strongly recommend to use at least Boost 1.67.

FFTW
    For some algorithms (P\ :math:`^3`\ M), |es| needs the FFTW library
    version 3 or later  [5]_ for Fourier transforms, including header
    files.

MPI
    Because |es| is parallelized with MPI, you need a working MPI
    environment that implements the MPI standard version 1.2.

Python
    |es|'s main user interface is via the Python 3 scripting interface.

Cython
    Cython is used for connecting the C++ core to Python.
    At least version 0.23 is required.


.. _Installing requirements on Ubuntu Linux:

Installing requirements on Ubuntu Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To make |es| run on 18.04 LTS, its dependencies can be installed with:

.. code-block:: bash

    sudo apt install build-essential cmake cython3 python3-numpy \
      libboost-all-dev openmpi-common fftw3-dev libhdf5-dev libhdf5-openmpi-dev \
      python3-opengl libgsl-dev


Optionally the ccmake utility can be installed for easier configuration:

.. code-block:: bash

    sudo apt install cmake-curses-gui

To run the tutorials and generate the documentation, additional Python packages
are required:

.. code-block:: bash

    sudo apt install python3-matplotlib python3-scipy ipython3 jupyter-notebook
    sudo pip3 install 'pint>=0.9'

Nvidia GPU acceleration
"""""""""""""""""""""""

If your computer has an Nvidia graphics card, you should also download and install the
CUDA SDK to make use of GPU computation:

.. code-block:: bash

    sudo apt install nvidia-cuda-toolkit

On Ubuntu, the default GCC compiler is too recent for nvcc, which will generate
compiler errors. You can either install an older version of GCC and select it
with environment variables ``CC`` and ``CXX`` when building |es|, or edit the
system header files as shown in the following example for Ubuntu 18.04:

.. code-block:: bash

    sudo sed -i 's/__GNUC__ > 6/__GNUC__ > 7/g' /usr/include/crt/host_config.h
    sudo sed -i 's/than 6/than 7/g' /usr/include/crt/host_config.h

AMD GPU acceleration
""""""""""""""""""""

If your computer has an AMD graphics card, you should also download and install the
ROCm SDK to make use of GPU computation:

.. code-block:: bash

    wget -qO - http://repo.radeon.com/rocm/apt/debian/rocm.gpg.key | sudo apt-key add -
    echo 'deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ xenial main' | sudo tee /etc/apt/sources.list.d/rocm.list
    sudo apt update
    sudo apt install libnuma-dev rocm-dkms rocblas rocfft rocrand rocthrust

After installing the ROCm SDK, please reboot your computer.


.. _Installing requirements on other Linux distributions:

Installing requirements on other Linux distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the following Dockerfiles to find the minimum set of packages
required to compile |es| on other Linux distributions:

* `CentOS <https://github.com/espressomd/docker/blob/master/docker/Dockerfile-centos>`_
* `Fedora <https://github.com/espressomd/docker/blob/master/docker/Dockerfile-fedora>`_
* `Debian <https://github.com/espressomd/docker/blob/master/docker/Dockerfile-debian>`_
* `OpenSUSE <https://github.com/espressomd/docker/blob/master/docker/Dockerfile-opensuse>`_


.. _Installing requirements on Mac OS X:

Installing requirements on Mac OS X
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Preparation
"""""""""""

To make |es| run on Mac OS X 10.9 or higher, you need to install its
dependencies. There are two possibilities for this, MacPorts and Homebrew.
We recommend MacPorts, but if you already have Homebrew installed, you can use
that too. To check whether you already have one or the other installed, run the
following commands:

.. code-block:: bash

    test -e /opt/local/bin/port && echo "MacPorts is installed"
    test -e /usr/local/bin/brew && echo "Homebrew is installed"

If both are installed, you need to remove one of the two. To do that, run one
of the following two commands:

.. code-block:: bash

    sudo port -f uninstall installed && rm -r /opt/local
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"

If Homebrew is already installed, you should resolve any problems reported by
the command

.. code-block:: bash

    brew doctor

If Anaconda Python or the Python from www.python.org are installed, you
will likely not be able to run |es|. Therefore, please uninstall them
using the following commands:

.. code-block:: bash

    sudo rm -r ~/anaconda[23]
    sudo rm -r /Library/Python

If you want to install MacPorts, download the installer package
appropriate for your Mac OS X version from
https://www.macports.org/install.php and install it.

If you want to install Homebrew, use the following commands.

.. code-block:: bash

    sudo xcode-select --install
    sudo xcodebuild -license accept
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Installing packages using MacPorts
""""""""""""""""""""""""""""""""""

Run the following commands:

.. code-block:: bash

    sudo xcode-select --install
    sudo xcodebuild -license accept
    sudo port selfupdate
    sudo port install cmake python37 py37-cython py37-numpy \
      openmpi-default fftw-3 +openmpi boost +openmpi +python37 \
      doxygen py37-opengl py37-sphinx gsl hdf5 +openmpi \
      py37-matplotlib py37-ipython py37-jupyter
    sudo port select --set cython cython37
    sudo port select --set python3 python37
    sudo port select --set mpi openmpi-mp


Installing packages using Homebrew
""""""""""""""""""""""""""""""""""

.. code-block:: bash

    brew install cmake python cython boost boost-mpi fftw \
      doxygen gsl numpy ipython jupyter
    brew install hdf5
    brew link --force cython
    pip install PyOpenGL matplotlib

.. _Quick installation:

Quick installation
------------------

If you have installed the requirements (see section :ref:`Requirements`) in
standard locations, compiling |es| is usually only a matter of creating a build
directory and calling ``cmake`` and ``make`` in it. See for example the command
lines below (optional steps which modify the build process are commented out):

.. code-block:: bash

    mkdir build
    cd build
    #cp myconfig-default.hpp myconfig.hpp # use the default configuration as template
    #nano myconfig.hpp                    # edit to add/remove features as desired
    cmake ..
    #ccmake . // in order to add/remove features like ScaFaCoS or CUDA
    make

This will build |es| with a default feature set, namely
:file:`src/config/myconfig-default.hpp`. This file is a C++ header file,
which defines the features that should be compiled in.
You may want to adjust the feature set to your needs. This can be easily done
by copying the :file:`myconfig-sample.hpp` which has been created in the :file:`build`
directory to :file:`myconfig.hpp` and only uncomment the features you want to use in your simulation.

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
command:

.. code-block:: bash

    ./pypresso <SCRIPT>

where ``<SCRIPT>`` is a ``python`` script which has to
be written by the user. You can find some examples in the :file:`samples`
folder of the source code directory. If you want to run in parallel, you should
have compiled with *Open MPI*, and need to tell MPI to run in parallel. The actual
invocation is implementation dependent, but in many cases, such as
*Open MPI*, you can use

.. code-block:: bash

    mpirun -n <N> ./pypresso <SCRIPT>

where ``<N>`` is the number of processors to be used.


.. _Features:

Features
--------

This chapter describes the features that can be activated in |es|. Even if
possible, it is not recommended to activate all features, because this
will negatively effect |es|'s performance.

Features can be activated in the configuration header :file:`myconfig.hpp` (see
section :ref:`myconfig.hpp\: Activating and deactivating features`). To
activate ``FEATURE``, add the following line to the header file:

.. code-block:: c++

    #define FEATURE


.. _General features:

General features
^^^^^^^^^^^^^^^^

-  ``ELECTROSTATICS`` This enables the use of the various electrostatics algorithms, such as P3M.

   .. seealso:: :ref:`Electrostatics`

-  ``MMM1D_GPU``

-  ``_P3M_GPU_FLOAT``


-  ``DIPOLES`` This activates the dipole-moment property of particles; In addition,
   the various magnetostatics algorithms, such as P3M are switched on.

   .. seealso::

       :ref:`Magnetostatics / Dipolar interactions`
       :ref:`Electrostatics`

-  ``SCAFACOS_DIPOLES``

-  ``ROTATION`` Switch on rotational degrees of freedom for the particles, as well as
   the corresponding quaternion integrator.

   .. seealso:: :ref:`Setting up particles`

   .. note::
      When this feature is activated, every particle has three
      additional degrees of freedom, which for example means that the
      kinetic energy changes at constant temperature is twice as large.

-  ``LANGEVIN_PER_PARTICLE`` Allows to choose the Langevin temperature and friction coefficient
   per particle.

-  ``ROTATIONAL_INERTIA``

-  ``EXTERNAL_FORCES`` Allows to define an arbitrary constant force for each particle
   individually. Also allows to fix individual coordinates of particles,
   keep them at a fixed position or within a plane.

-  ``MASS`` Allows particles to have individual masses. Note that some analysis
   procedures have not yet been adapted to take the masses into account
   correctly.

   .. seealso:: :attr:`espressomd.particle_data.ParticleHandle.mass`

-  ``EXCLUSIONS`` Allows to exclude specific short ranged interactions within
   molecules.

   .. seealso:: :meth:`espressomd.particle_data.ParticleHandle.add_exclusion`

-  ``COMFIXED`` Allows to fix the center of mass of all particles of a certain type.

-  ``BOND_CONSTRAINT`` Turns on the RATTLE integrator which allows for fixed lengths bonds
   between particles.

-  ``VIRTUAL_SITES_RELATIVE`` Virtual sites are particles, the position and velocity of which is
   not obtained by integrating equations of motion. Rather, they are
   placed using the position (and orientation) of other particles. The
   feature allows for rigid arrangements of particles.

   .. seealso:: :ref:`Virtual sites`

-  ``COLLISION_DETECTION`` Allows particles to be bound on collision.

In addition, there are switches that enable additional features in the
integrator or thermostat:

-  ``NPT`` Enables an on-the-fly NPT integration scheme.

   .. seealso:: :ref:`Isotropic NPT thermostat`


-  ``REACTION_ENSEMBLE``

-  ``ENGINE``

-  ``PARTICLE_ANISOTROPY``


.. _Fluid dynamics and fluid structure interaction:

Fluid dynamics and fluid structure interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``DPD`` Enables the dissipative particle dynamics thermostat and interaction.

   .. seealso:: :ref:`DPD interaction`

-  ``LB_BOUNDARIES``

-  ``LB_BOUNDARIES_GPU``

-  ``LB_ELECTROHYDRODYNAMICS`` Enables the implicit calculation of electro-hydrodynamics for charged
   particles and salt ions in an electric field.

-  ``ELECTROKINETICS``

-  ``EK_BOUNDARIES``

-  ``EK_DEBUG``

-  ``EK_DOUBLE_PREC``


.. _Interaction features:

Interaction features
^^^^^^^^^^^^^^^^^^^^

The following switches turn on various short ranged interactions (see
section :ref:`Isotropic non-bonded interactions`):

-  ``TABULATED`` Enable support for user-defined non-bonded interaction potentials.

-  ``LENNARD_JONES`` Enable the Lennard-Jones potential.

-  ``LENNARD_JONES_GENERIC`` Enable the generic Lennard-Jones potential with configurable
   exponents and individual prefactors for the two terms.

-  ``LJCOS`` Enable the Lennard-Jones potential with a cosine-tail.

-  ``LJCOS2`` Same as ``LJCOS``, but using a slightly different way of smoothing the
   connection to 0.

-  ``GAY_BERNE`` (experimental)

-  ``HERTZIAN``

-  ``NO_INTRA_NB``

-  ``MORSE`` Enable the Morse potential.

-  ``BUCKINGHAM`` Enable the Buckingham potential.

-  ``SOFT_SPHERE`` Enable the soft sphere potential.

-  ``SMOOTH_STEP`` Enable the smooth step potential, a step potential with two length
   scales.

-  ``BMHTF_NACL`` Enable the Born-Meyer-Huggins-Tosi-Fumi potential, which can be used
   to model salt melts.

-  ``GAUSSIAN``

-  ``HAT``

-  ``UMBRELLA`` (experimental)

Some of the short-range interactions have additional features:

-  ``LJGEN_SOFTCORE`` This modifies the generic Lennard-Jones potential
   (``LENNARD_JONES_GENERIC``) with tunable parameters.


.. _Debug messages:

Debug messages
^^^^^^^^^^^^^^

Finally, there is a flag for debugging:

-  ``ADDITIONAL_CHECKS`` Enables numerous additional checks which can detect
   inconsistencies especially in the cell systems. These checks are however
   too slow to be enabled in production runs.

   .. note::
      Because of a bug in OpenMPI versions 2.0-2.1, 3.0.0-3.0.2 and 3.1.0-3.1.2
      that causes a segmentation fault when running the |es| OpenGL visualizer
      with feature ``ADDITIONAL_CHECKS`` enabled together with either
      ``ELECTROSTATICS`` or ``DIPOLES``, the subset of additional checks for
      those two features are disabled if an unpatched version of OpenMPI is
      detected during compilation.


Features marked as experimental
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some of the above features are marked as EXPERIMENTAL. Activating these features can have unexpected side effects and some of them have known issues. If you activate any of these features, you should understand the corresponding source code and do extensive testing. Furthermore, it is necessary to define ``EXPERIMENTAL_FEATURES`` in :file:`myconfig.hpp`.


External features
^^^^^^^^^^^^^^^^^

External features cannot be added to the :file:`myconfig.hpp` file by the user.
They are added by CMake if the corresponding dependency was found on the
system. Some of these external features are optional and must be activated
using a CMake flag (see :ref:`Options and Variables`).

- ``CUDA`` Enables GPU-specific features.

- ``FFTW`` Enables features relying on the fast Fourier transforms, e.g. P3M.

- ``H5MD`` Write data to H5MD-formatted hdf5 files (see :ref:`Writing H5MD-files`)

- ``SCAFACOS`` Enables features relying on the ScaFaCoS library (see
  :ref:`ScaFaCoS electrostatics`, :ref:`ScaFaCoS magnetostatics`).

- ``GSL`` Enables features relying on the GNU Scientific Library, e.g.
  :meth:`espressomd.cluster_analysis.Cluster.fractal_dimension`.



.. _Configuring:

Configuring
-----------

.. _myconfig.hpp\: Activating and deactivating features:

:file:`myconfig.hpp`: Activating and deactivating features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|es| has a large number of features that can be compiled into the binary.
However, it is not recommended to actually compile in all possible
features, as this will slow down |es| significantly. Instead, compile in only
the features that are actually required. A strong gain in speed can be
achieved by disabling all non-bonded interactions except for a single
one, e.g. ``LENNARD_JONES``. For developers, it is also possible to turn on or off a
number of debugging messages. The features and debug messages can be
controlled via a configuration header file that contains C-preprocessor
declarations. Subsection :ref:`Features` describes all available features. If a
file named :file:`myconfig.hpp` is present in the build directory when ``cmake``
is run, all features defined in it will be compiled in. If no such file exists,
the configuration file :file:`src/config/myconfig-default.hpp` will be used
instead, which turns on the default features.

When you distinguish between the build and the source directory, the
configuration header can be put in either of these. Note, however, that
when a configuration header is found in both directories, the one in the
build directory will be used.

By default, the configuration header is called :file:`myconfig.hpp`.
The configuration header can be used to compile different binary
versions of with a different set of features from the same source
directory. Suppose that you have a source directory :file:`$srcdir` and two
build directories :file:`$builddir1` and :file:`$builddir2` that contain
different configuration headers:

*  :file:`$builddir1/myconfig.hpp`:

  .. code-block:: c++

    #define ELECTROSTATICS
    #define LENNARD_JONES

*  :file:`$builddir2/myconfig.hpp`:

  .. code-block:: c++

    #define LJCOS

Then you can simply compile two different versions of |es| via:

.. code-block:: bash

    cd builddir1
    cmake ..
    make

    cd builddir2
    cmake ..
    make

To see what features were activated in :file:`myconfig.hpp`, run:

.. code-block:: bash

    ./pypresso

and then in the Python interpreter:

.. code-block:: python

    import espressomd
    print(espressomd.features())


.. _cmake:

``cmake``
^^^^^^^^^

In order to build the first step is to create a build directory in which
cmake can be executed. In cmake, the *source directory* (that contains
all the source files) is completely separated from the *build directory*
(where the files created by the build process are put). ``cmake`` is
designed to *not* be executed in the source directory. ``cmake`` will
determine how to use and where to find the compiler, as well as the
different libraries and tools required by the compilation process. By
having multiple build directories you can build several variants of |es|,
each variant having different activated features, and for as many
platforms as you want.

Once you've run ``ccmake``, you can list the configured variables with
``cmake -LAH -N | less`` (uses a pager) or with ``ccmake ..`` and pressing
key ``t`` to toggle the advanced mode on (uses the curses interface).

**Example:**

When the source directory is :file:`srcdir` (the files where unpacked to this
directory), then the user can create a build directory :file:`build` below that
path by calling :file:`mkdir srcdir/build`. In the build directory ``cmake`` is to be
executed, followed by a call to make. None of the files in the source directory
are ever modified by the build process.

.. code-block:: bash

    cd build
    cmake ..
    make

Afterwards |es| can be run via calling :file:`./pypresso` from the command line.


.. _ccmake:

``ccmake``
^^^^^^^^^^

Optionally and for easier use, the curses interface to cmake can be used
to configure |es| interactively.

**Example:**

Alternatively to the previous example, instead of cmake, the ccmake executable
is called in the build directory to configure |es|, followed by a call to make:

.. code-block:: bash

    cd build
    ccmake ..
    make

Fig. :ref:`ccmake-figure` shows the interactive ccmake UI.

.. _ccmake-figure:

.. figure:: figures/ccmake-example.png
   :alt: ccmake interface
   :width: 70.0%
   :align: center

   ccmake interface


.. _Options and Variables:

Options and Variables
^^^^^^^^^^^^^^^^^^^^^

The behavior of |es| can be controlled by means of options and variables
in the :file:`CMakeLists.txt` file. Also options are defined there. The following
options are available:

* ``WITH_CUDA``: Build with GPU support

* ``WITH_HDF5``: Build with HDF5

* ``WITH_TESTS``: Enable tests

* ``WITH_SCAFACOS``: Build with ScaFaCoS support

* ``WITH_VALGRIND_INSTRUMENTATION``: Build with valgrind instrumentation
  markers

When the value in the :file:`CMakeLists.txt` file is set to ON, the corresponding
option is created; if the value of the option is set to OFF, the
corresponding option is not created. These options can also be modified
by calling ``cmake`` with the command line argument ``-D``:

.. code-block:: bash

    cmake -D WITH_HDF5=OFF srcdir

When an option is activated, additional options may become available.
For example with ``-D WITH_CUDA=ON``, one can choose the CUDA compiler with
``-D WITH_CUDA_COMPILER=<compiler_id>``, where ``<compiler_id>`` can be
``nvcc`` (default), ``clang`` or ``hip``. For ``hip``, an additional
``-D ROCM_HOME=<path_to_rocm>`` variable becomes available, with default value
``ROCM_HOME=/opt/rocm``.

Environment variables can be passed to CMake. For example, to select Clang, use
``CC=clang CXX=clang++ cmake .. -DWITH_CUDA=ON -DWITH_CUDA_COMPILER=clang``.
If you have multiple versions of the CUDA library installed, you can select the
correct one with ``CUDA_BIN_PATH=/usr/local/cuda-10.0 cmake .. -DWITH_CUDA=ON``
(with Clang as the CUDA compiler, you also need to override its default CUDA
path with ``-DCMAKE_CXX_FLAGS=--cuda-path=/usr/local/cuda-10.0``).


Compiling, testing and installing
---------------------------------

The command ``make`` is mainly used to compile the source code, but it
can do a number of other things. The generic syntax of the ``make``
command is:

.. code-block:: bash

    make [options] [target] [variable=value]

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
    Install |es| in the path specified by the CMake variable
    ``CMAKE_INSTALL_PREFIX``. The path can be changed by calling CMake
    with ``cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/espresso``. Do not use
    ``make DESTDIR=/path/to/espresso install`` to install to a specific path,
    this will cause issues with the runtime path (RPATH) and will conflict
    with the CMake variable ``CMAKE_INSTALL_PREFIX`` if it has been set.

``doxygen``
    Creates the Doxygen code documentation in the :file:`doc/doxygen`
    subdirectory.

``sphinx``
    Creates the ``sphinx`` code documentation in the :file:`doc/sphinx`
    subdirectory.

``tutorials``
    Creates the tutorials in the :file:`doc/tutorials` subdirectory.

``doc``
    Creates all documentation in the :file:`doc` subdirectory (only when
    using the development sources).

A number of options are available when calling ``make``. The most
interesting option is probably ``-j num_jobs``, which can be used for
parallel compilation on computers that have more than one CPU or core.
*num_jobs* specifies the maximal number of jobs that will be run.
Setting *num_jobs* to the number of available processors speeds up the
compilation process significantly.

.. _Running es:

Running |es|
------------

Executing a simulation script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|es| is implemented as a Python module. This means that you need to write a
python script for any task you want to perform with |es|. In this chapter,
the basic structure of the interface will be explained. For a practical
introduction, see the tutorials, which are also part of the
distribution. To use |es|, you need to import the espressomd module in your
Python script. To this end, the folder containing the python module
needs to be in the Python search path. The module is located in the
:file:`src/python` folder under the build directory. A convenient way to run
python with the correct path is to use the pypresso script located in
the build directory.

.. code-block:: bash

    ./pypresso simulation.py

The ``pypresso`` script is just a wrapper in order to expose the |es| python
module to the system's python interpreter by modifying the ``$PYTHONPATH``.
Please see the following chapter :ref:`Setting up the system` describing how
to actually write a simulation script for |es|.

Running an interactive notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running the Jupyter interpreter requires using the ``ipypresso`` script, which
is also located in the build directory (its name comes from the IPython
interpreter, today known as Jupyter). To run the tutorials, you will need
to start the Jupyter interpreter in notebook mode:

.. code-block:: bash

    cd doc/tutorials
    ../../ipypresso notebook

You may then browse through the different tutorial folders. Files whose name
ends with extension .ipynb can be opened in the browser. Click on the Run
button to execute the current block, or use the keyboard shortcut Shift+Enter.
If the current block is a code block, the ``In [ ]`` label to the left will
change to ``In [*]`` while the code is being executed, and become ``In [1]``
once the execution has completed. The number increments itself every time a
code cell is executed. This bookkeeping is extremely useful when modifying
previous code cells, as it shows which cells are out-of-date. It's also
possible to run all cells by clicking on the "Run" drop-down menu, then on
"Run All Below". This will change all labels to ``In [*]`` to show that the
first one is running, while the subsequent ones are awaiting execution.
You'll also see that many cells generate an output. When the output becomes
very long, Jupyter will automatically put it in a box with a vertical scrollbar.
The output may also contain static plots, dynamic plots and videos. It is also
possible to start a 3D visualizer in a new window, however closing the window
will exit the Python interpreter and Jupyter will notify you that the current
Python kernel stopped. If a cell takes too long to execute, you may interrupt
it with the stop button.

To close the Jupyter notebook, go to the terminal where it was started and use
the keyboard shortcut Ctrl+C twice.

When starting the Jupyter interpreter in notebook mode, you may see the
following warning in the terminal:

.. code-block:: none

    [TerminalIPythonApp] WARNING | Subcommand `ipython notebook` is deprecated and will be removed in future versions.
    [TerminalIPythonApp] WARNING | You likely want to use `jupyter notebook` in the future

This only means |es| was compiled with IPython instead of Jupyter. If Jupyter
is installed on your system, the notebook will automatically close IPython and
start Jupyter. To recompile |es| with Jupyter, provide ``cmake`` with the flag
``-DIPYTHON_EXECUTABLE=$(which jupyter)``.

You can find the official Jupyter documentation at
https://jupyter.readthedocs.io/en/latest/running.html

.. _Debugging es:

Debugging |es|
--------------

Exceptional situations occur in every program.  If |es| crashes with a
segmentation fault, that means that there was a memory fault in the
simulation core which requires running the program in a debugger.  The
``pypresso`` executable file is actually not a program but a script
which sets the Python path appropriately and starts the Python
interpreter with your arguments.  Thus it is not possible to directly
run ``pypresso`` in a debugger.  However, we provide some useful
command line options for the most common tools.

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
