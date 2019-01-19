
.. _Contact the Developers:

Contact the Developers
======================

To contact the |es| developers, please write an email to the developers mailing list:
espressomd-devel@nongnu.org
to subscribe to the developers' mailing list go to
http://lists.nongnu.org/mailman/listinfo/espressomd-devel


.. _Before you start a development project:

Before you start a development project
--------------------------------------
Before you start a development project for |es|, please always write to the developers mailing list and describe the project.
This is to avoid that several people work on the same thing at the same time. Also, implementation details can be discussed in advance. In many cases, existing developers can point to re-usable code and simpler solutions.


.. _Development Environment:

Development Environment
=======================


.. _Required Development Tools:

Required Development Tools
--------------------------

-  First of all, please install the dependencies for compiling |es|. See the section on "Getting, compiling and running" in the user guide.

-  To be able to access the development version of |es|, you will need
   the distributed versioning control system git_.

-  To build the sphinx documentation, you will need the Python packages listed in :file:`requirements.txt` in the top-level source directory. To install them, issue:

   .. code-block:: bash

      pip install --upgrade --user -r requirements.txt

   Note, that some distributions now use ``pip`` for Python3 and ``pip2`` for Python 2.

-  To build the tutorials, you will need LaTeX.

-  To compile the Doxygen code documentation, you will need to have the
   tool doxygen_.

All of these tools should be easy to install on most Unix operating
systems.

.. _Getting the Development Code:

Getting the Development Code
----------------------------
We use Github for storing the source code and its history, and for managing the development process.
The repository is located at http://github.com/espressomd/espresso.
To get the current development code, run:

.. code-block:: bash

  git clone git://github.com/espressomd/espresso

This will create a directory named "espresso" which contains the code.
The build process does not differ from the one for release versions described in the users' guide.


Build System
------------

The build system of |es| is based on CMake.

The central source files of the build system are the following:

-  :file:`CMakeLists.txt`

-  Contents of the :file:`cmake` directory

-  The :file:`CMakeLists.txt` files in the :file:`src/`, :file:`doc/`, and :file:`testsuite/` directories and their sub-directories

The most common reasons for editing these files are:

-  Adding new source files
-  Adding new external dependencies

Adding New Source Files
~~~~~~~~~~~~~~~~~~~~~~~

To add new files to |es| (like C++ source files or header files) you
need to look at the CMakeLists.txt in the directory where the file is located.

* Please note that .hpp-header files usually do not have to be added to CMakeList.txt

* All files are explicitly included (e.g., testsuite/python/CMakeLists.txt):: 

      set(py_tests  bondedInteractions.py
                   cellsystem.py
                   constraint_shape_based.py
                   coulomb_cloud_wall.py)


Testsuite
---------

-  New or significantly changed features will only be accepted if they have a test case.
   This is to make sure the feature is not broken by future changes to |es|, and so other users can get an impression of what behavior is guaranteed to work.
-  There are multiple kinds of tests:

  -  C++-unit tests, testing individual C++ functions and classes. They make use of the boost unit test framework and reside in :file:`/src/core/unit_tests`
  -  Python integration tests, testing the Python interface and (physical) results of features. They reside in :file:`/testsuite/python`
  -  CMake tests, testing the software can be successfully installed. They reside in :file:`/testsuite/cmake`
  -  Python scripts tests, testing the IPython notebooks and Python samples in :file:`/doc/tutorials` and :file:`/samples`. They reside in :file:`/testsuite/scripts/tutorials` resp. :file:`/testsuite/scripts/samples` and are executed on a different schedule

- To execute the tests, run:

  .. code-block:: bash

     make check

  in the top build directory.

- See :ref:`Unit testing` for how to develop new tests


.. _Documentation:

Documentation
=============

The documentation of |es| consists of four parts:

  -  The users' guide and developers' guide are located in :file:`doc/sphinx`, and make use of the Sphinx Python package
  -  In-code documentation for the Python interface is located in the various files in src/python/espressomd and also makes use of the Sphinx Python package. We make use of the napoleon extension and use the NumPy documentation style.
  -  In-code documentation of the C++ core is located in the .cpp and .hpp files in :file:`/src/core` and its sub-directories and makes use of Doxygen.

Doxygen Code Documentation
--------------------------

The documentation of each function should contain a short description,
if necessary a more detailed description and a description for the
return value and parameters.

Look at the documentation of existing files and functions to get a
feeling how it should be!

Doxygen is able to understand simple LaTeX and HTML commands as well as
some special command in order to give the documentation a nice structure
and to make it more readable. In the following list you find a short
description of the most common commands we need:

-  | ``\anchor`` *name* *description*
   | Create an anchor to which you can refer using the ``\ref`` command.

-  | ``\ref`` *name* ``["``\ *text*\ ``"]``
   | Insert a link to another object in the documentation (*e.g.*\ an
     anchor).

-  | ``<a href="http://www.your_url.html">title</a>``
   | Link to an external HTML source.

-  | ``\file`` *name* *description*
   | Special anchor for a file.

-  | ``\image html`` *image*
   | Include a picture. The picture file should reside in the subdir
     ``doc/doxygen/figs``. Do not use the HTML ``<img>``-tag to include
     pictures, as doxygen_ will not copy the pictures into the
     documentation.

-  | ``<ul> <li>List entry 1</li> <li>List entry 2</li></ul>``
   | Creates a list in the documentation.

-  | ``\param`` *name* *description*
   | Document the parameter of a function.

-  | ``\return`` *description*
   | Document the return value of a function.

.. _Programmers's Guide:


Programmer's Guide
==================

This chapter provides some hints on how to extend |es|. It is not
exhaustive, so for major changes the best documentation are the other
developers.


Source code structure
---------------------

The source tree has the following structure:

* src: The actual source code

  * core: The C++ source code of the simulation core
  * python/espressomd: Source of the espressomd Python module and its submodules
  * script_interface: C++ source code of the script_interface component, which links Python classes to functionality in the simulation core

* doc: Documentation

  * sphinx: The sphinx-based documentation, consisting of user and developer guide.
  * tutorials/python: Source and pdf files for the introductory tutorials
  * doxygen_: Build directory for the C++ in-code documentation

* testsuite/python: Python integration tests. Note that some C++ unit tests for individual core components are in src/core/unittests
* samples/python: Some sample scripts
* libs: External dependencies (at this point h5xx)
* maintainer: Files used by the maintainers

  * configs: Collection of myconfig.hpp files which activate different sets of features for testing.
  * docker: Definitions of the docker images for various distributions used for continuous integration testing
  * CI: Support files for the continuous integration testing run on the Travis-CI service.
  * jenkins: Outdated support files for the Jenkins continuous integration testing


Flow control and communications architecture
--------------------------------------------
Espresso uses two communication models, namely master-slave and synchronous.

* When Espresso does not run an integration, it works in the master-slave mode, i.e. the head node (MPI rank 0) in a parallel simulation
  runs the Python script, whereas all other nodes are idle until they receive a command from the head node. Such commands include particle creation,
  changing of particle properties and changing global simulation parameters.
  When a Python command such as::

    system.part.add(pos=(1, 2, 3))

  is issued, the head node determines, which node is responsible for the given position, and then sends the node the command to place the particle.

* When an integration is started in Python on the head node, a command to start the integration is sent to all nodes, in the master-slave framework described above.
  Then, Espresso switches into the synchronous mode, in which all nodes run the same code in the integration loop at the same time.
  The code of the main integration loop is in ``integrate.cpp:integrate_vv()``.
  When writing code which is run during the main integration loop, no commands making use of the master-slave mechanism can be called.
  When code during the integration loop executes MPI communication, it has to be ensured, that the MPI call is executed on all nodes
  involved in the communication. If this is not done, a deadlock will result.

Adding calls to the master-slave framework
------------------------------------------

Using an instance of MpiCallback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Write the callback slave function, which will be executed on all nodes except the head node (0):

  .. code-block:: c++

    void my_callback(int p1, int p2) {
      // Do something. The two int-parameters can be used for anything
    }

* On all nodes, the callback has to be registered:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add(my_callback);
    }

  You can, e.g., call your registration from ``initialize.cpp:on_program_start()``
  Instead of a static function, from which a ``std::function<void(int,int)>`` can be constructed can
  be used. For example:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add([](int, int){ /* Do something */ });
    }

  can be used to add a lambda function as callback.
* Then, you can use your callback from the head node:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void call_my_callback() {
      Communication::mpiCallbacks.call(my_callback, param1, param2);
    }

  This only works outside the integration loop. After the callback has been called, synchronous mpi communication can be done.

Legacy callbacks
~~~~~~~~~~~~~~~~

Older code uses callbacks defined in the ``CALLBACK_LIST`` preprocessor macro in :file:`communications.cpp`. They are called via ``mpi_call()``.
See ``communications.cpp:mpi_place_particle()`` for an example.

Adding New Bonded Interactions
------------------------------

To add a new bonded interaction, the following steps have to be taken

* Simulation core:

  * Define a structure holding the parameters (prefactors, etc.) of the interaction
  * Write functions for calculating force and energy, respectively.
  * Write a setter function, which takes the parameters of the interactions and stores them in the bonded interactions data structure
  * Add calls to the force and energy calculation functions to the force calculation in the integration loop as well as to energy and pressure/stress tensor analysis

* Python interface

  * Import the definition of the bond data structure from the simulation core
  * Implement a class for the bonded interaction derived from the BondedInteraction base class

Defining the data structure for the interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The data structures for bonded interactions reside in ``src/core/bonded_interactions/bonded_interaction_data.hpp``.

* Add your interaction to the ``enum BondedInteraction``.
  This enumeration is used to identify different bonded interactions.
* Add a typedef struct containing the parameters of the interaction. Use the one for the FENE interaction as template:

  .. code-block:: c++

    typedef struct {
      double k;
      [...]
    } Fene_bond_parameters;

* Add a member to the typedef union Bond_parameters. For the FENE bond it looks like this:

  .. code-block:: c++

    Fene_bond_parameters fene;


Functions for calculating force and energy, and for setting parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every interaction resides in its own source .cpp and .hpp. A simple example for a
bonded interaction is the FENE bond in :file:`src/core/bonded_interactions/fene.cpp` and :file:`src/core/bonded_interactions/fene.hpp`.
Use these two files as templates for your interaction.

Notes:

* The names of function arguments mentioned below are taken from the FENE bond in :file:`src/core/bonded_interactions/fene.cpp` and :file:`src/core/bonded_interactions/fene.hpp`. It is recommended to use the same names for the corresponding functions for your interaction.
* The recommended signatures of the force and energy functions are:

  .. code-block:: c++

    inline int calc_fene_pair_force(Particle *p1, Particle *p2,
                                    Bonded_ia_parameters *iaparams,
                                    double dx[3], double force[3])
    inline int fene_pair_energy(Particle *p1, Particle *p2,
                                Bonded_ia_parameters *iaparams,
                                double dx[3], double *_energy)

  Here, ``fene`` needs to be replaced by the name of the new interaction.
* The setter function gets a ``bond_type`` which is a numerical id identifying the number of the bond type in the simulation. It DOES NOT determine the type of the bond potential (harmonic vs FENE).
  The signature of the setter function has to contain the ``bond_type``, the remaining parameters are specific to the interaction. For the FENE bond, e.g., we have:

  .. code-block:: c++

    fene_set_params(int bond_type, double k, double drmax, double r0)

  A return value of ``ES_OK`` is returned on success, ``ES_ERR`` on error, e.g., when parameters are invalid.
* The setter function must call ``make_bond_type_exists()`` with that bond type, to allocate the memory for storing the parameters.
* Afterwards, the bond parameters can be stored in the global variable ``bonded_ia_params[bond_type]``

  * ``bonded_ia_params[bond_type].num`` is the number of particles involved in the bond -1. I.e., 1 for a pairwise bonded potential such as the FENE bond.
  * The parameters for the individual bonded interaction go to the member of ``Bond_parameters`` for your interaction defined in the previous step. For the FENE bond, this would be:

    .. code-block:: c++

      bonded_ia_params[bond_tpe].p.fene

* At the end of the parameter setter function, do not forget the call to ``mpi_bcast_ia_params()``, which will sync the parameters just set to other compute nodes in a parallel simulation.
* The routines for calculating force and energy return an integer. A return value of 0 means OK, a value of 1 means that the particles are too far apart and the bond is broken. This will stop the integration with a runtime error.
* The functions for calculating force and energy can make use of a pre-calculated distance vector (dx) pointing from particle 2 to particle 1.
* The force on particle 1 has to be stored in the force vector  (not added to it). The force on particle 2 will be obtained from Newton's law.
* The result of the energy calculation is placed in (NOT added to) the ``_energy`` argument of the energy calculation function.



Including the bonded interaction in the force calculation and the energy and pressure analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* In :file:`src/core/bonded_interactions/bonded_interaction_data.cpp`:

    #. Add a name for the interaction to ``get_name_of_bonded_ia()``.
    #. In ``calc_maximal_cutoff()``, add a case for the new interaction which
       makes sure that ``max_cut`` is larger than the interaction range of the
       new interaction, typically the bond length.  This is necessary to ensure
       that, in a parallel simulation, a compute node has access to both bond
       partners. This value is always used as calculated by
       ``calc_maximal_cutoff``, therefore it is not strictly necessary that the
       maximal interaction range is stored explicitly.
    #. Besides this, you have enter the force respectively the energy
       calculation routines in ``add_bonded_force``, ``add_bonded_energy``,
       ``add_bonded_virials`` and ``pressure_calc``. The pressure occurs ice,
       once for the parallelized isotropic pressure and once for the tensorial
       pressure calculation. For pair forces, the pressure is calculated using
       the virials, for many body interactions currently no pressure is
       calculated.
    #. Do not forget to include the header file of your interaction.

* Force calculation: in :file:`forces_inline.hpp` in the function
  ``add_bonded_force()``, add your bond to the switch statement. For the FENE
  bond, e.g., the code looks like this:

  .. code-block:: c++

    case BONDED_IA_FENE:
      bond_broken = calc_fene_pair_force(p1, p2, iaparams, dx, force);

* Energy calculation: add similar code to ``add_bonded_energy()`` in :file:`energy_inline.hpp`
* Pressure, stress tensor and virial calculation: If your bonded interaction is
  a pair bond and does not modify the particles involved, add similar code as
  above to ``pressure.hpp:calc_bonded_pair_force()``. Otherwise, you have to
  implement a custom solution for virial calculation.


Adding the bonded interaction in the Python interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please note that the following is Cython code (www.cython.org), rather than pure Python.

* In :file:`src/python/espressomd/interactions.pxd`:

  * import the parameter data structure from the C++ header file for your interaction. For the FENE bond, this looks like:

    .. code-block:: cython

      cdef extern from "interaction_data.hpp":
          ctypedef struct Fene_bond_parameters:
              double k
              double drmax
              double r0
              double drmax2
              double drmax2i

  * Add your bonded interaction to the Cython copy of the BondedInteractions enum analogous to the one in the core:, described above:

    .. code-block:: cython

      cdef enum enum_bonded_interaction "BondedInteraction":
          BONDED_IA_NONE = -1,
          BONDED_IA_FENE,
          BONDED_IA_HARMONIC,
          [...]

    The spelling has to match the one in the C++ enum exactly.
  * Adapt the Cython copy of the bond_parameters union analogous to the C++ core.  The member name has to match the one in C++ exactly:

    .. code-block:: cython

      ctypedef union bond_parameters "Bond_parameters":
          Fene_bond_parameters fene
          Oif_global_forces_bond_parameters oif_global_forces
          Oif_local_forces_bond_parameters oif_local_forces
          Harmonic_bond_parameters harmonic

  * Import the declaration of the setter function implemented in the core. For the FENE bond, this looks like:

    .. code-block:: cython

        cdef extern from "fene.hpp":
            int fene_set_params(int bond_type, double k, double drmax, double r0)

* In :file:`src/python/espressomd/interactions.pyx`:

  * Implement the Cython class for the bonded interaction, using the one for
    the FENE bond as template. Please use pep8 naming convention:

    .. code-block:: cython

        class FeneBond(BondedInteraction):

            def __init__(self, *args, **kwargs):
                """
                FeneBond initializer. Used to instantiate a FeneBond identifier
                with a given set of parameters.

                Parameters
                ----------
                k : float
                    Specifies the magnitude of the bond interaction.
                d_r_max : float
                          Specifies the maximum stretch and compression length of the
                          bond.
                r_0 : float, optional
                      Specifies the equilibrium length of the bond.
                """
                super(FeneBond, self).__init__(*args, **kwargs)

            def type_number(self):
                return BONDED_IA_FENE

            def type_name(self):
                return "FENE"

            def valid_keys(self):
                return "k", "d_r_max", "r_0"

            def required_keys(self):
                return "k", "d_r_max"

            def set_default_params(self):
                self._params = {"r_0": 0.}

            def _get_params_from_es_core(self):
                return \
                    {"k": bonded_ia_params[self._bond_id].p.fene.k,
                     "d_r_max": bonded_ia_params[self._bond_id].p.fene.drmax,
                     "r_0": bonded_ia_params[self._bond_id].p.fene.r0}

            def _set_params_in_es_core(self):
                fene_set_params(
                    self._bond_id, self._params["k"], self._params["d_r_max"], self._params["r_0"])

* In :file:`testsuite/python/bondedInteractions.py`:

  * Add a test case, which verifies that parameters set and gotten from the interaction are consistent::

        test_fene = generateTestForBondParams(
            0, FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})

.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
