=================
Developer's Guide
=================
.. warning::
   The information found in this version of the Developer's Guide is
   outdated.  Please see the section :ref:`Contact the Developers` and
   ask for advice if you plan to start working on |es|.


.. _Contact the Developers:

Contact the Developers
======================

To contact the |es| developers, please write an email to the developers mailing list:
espressomd-devel@nongnu.org
to subscribe to the developers’ mailing list go to
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

-  To build the sphinx documentation, you will need the Python packages listed in ``requirements.txt`` in the top-level source directory. To install them, issue::
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
The repository is located at
http://github.com/espressomd/espresso
To get the current development code, run
git clone git://github.com/espressomd/espresso
This will create a directory named "espresso" which contains the code.
The build process does not differ from the one for release versions described in the users' guide.


Build System
------------

The build system of |es| is based on CMake.

The central source files of the build system are the following:

-  ``CMakeList.txt``

-  Contents of the ``cmake`` directory

-  The CMakeList.txt files in the ``src/``, ``doc/``, and ``testsuite/`` directories and their sub-directories

The most common reasons for editing these files are:
-  Adding new source files
-  Adding new external dependencies

Adding New Source Files
~~~~~~~~~~~~~~~~~~~~~~~

To add new files to |es| (like C++ source files or header files) you
need to look at the CMakeList.txt in the directory where the file is located.
* Please note that .hpp-header files usually do not have to be added to CMakeList.txt
* In some cases (e.g., src/core/CMakeList.txt), the CMakeList.txt contains a wild-card include like this::

      file(GLOB EspressoCore_SRC *.cpp)

  In this case, placing a file with that ending is enough.

* In other cases, the files are explicitly included (e.g., testsuite/python/CMakeList):: 

      set(py_tests  bondedInteractions.py
                   cellsystem.py
                   constraint_shape_based.py
                   coulomb_cloud_wall.py)

  In that case, add the new file to the list.
   


Testsuite
---------

-  New or significantly changed features will only be accepted, if they have a test case. 
   This is to make sure, the feature is not broken by future changes to |es|, and so other users can get an impression of what behaviour is guaranteed to work.
-  There are two kinds of tests:

  -  C++-unit tests, testing individual C++ functions and classes. They make use of the boost unit test framework and reside in ``src/core/unit_tests``
  -  Python integration tests, testing the Python interface and (physical) results of features. They reside in ``testsuite/python``

-  To execute the tests, run::

     make check 

   in the top build directory.


.. _Documentation:

Documentation
=============

The documentation of |es| consists of four parts:

  -  The users' guide and developers' guide are located in ``doc/sphinx``, and make use of the Sphinx Python package
  -  In-code documentation for the Python interface is located in the various files in src/python/espressomd and also makes use of the Sphinx Python package. We make use of the napolean extension and use the NumPy documentation style.
  -  In-code documentation of the C++ core is located in the .cpp and .hpp files in ``/src/core`` and its sub-directories and makes use of Doxygen.



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

-  | ``\return`` *decription*
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
  * travis: Support files for the continuous integration testing run on the Travis-CI service.
  * jenkins: Outdated support files for the Jenkins continuous integration testing
		

Flow control and communications architecture
--------------------------------------------
Espresso uses two communication models, namely master-slave and synchronous.

* When Espresso does not run an integration, it works in the master-slave mode, i.e. the head node (MPI rank 0) in a parallel simulation
  runs the Python script, whereas all other nodes are idle until they receive a command from the head node. Such commands include particle creation,
  changing of particle properties and changing global simulation parameters.
  When a Python command such as:::

    system.part.add(pos=(1,2,3))

  is issued, the head node determines, which node is responsible for the given position, and then sends the node the command to place the particle.

* When an integration is started in Python on the head node, a command to start the integration is sent to all nodes, in the master-slave framework described above.
  Then, Espresso switches into the synchronous mode, in which all nodes run the same code in the integration loop at the same time.
  The code of the main integration loop is in integrate.cpp:integrate_vv().
  When writing code which is run during the main integration loop, no commands making use of the master-slave mechanism can be called.
  When code during the integration loop executes MPI communication, it has to be ensured, that the MPI call is executed on all nodes
  involved in the communication. If this is not done, a deadlock will result.

Adding calls to the master-slave framework
------------------------------------------

Using an instance of MpiCallback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Write the callback slave function, which will be executed on all nodes except the head node (0)::

    void my_callback(int p1, int p2) {
      // Do something. The two int-parameters can be usued for anything
    }

* On all nodes, the callback has to be registered::

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add(my_callback);
    }

  You can, e.g., call your registration from initialize.cpp:on_program_start()
  Instead of a static function, from which a ``std::function<void(int,int)>`` can be constructed can
  be used. For example::

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add([](int, int){ /* Do something */ });
    }

  can be used to add a lambda function as callback.
* Then, you can use your callback from the head node::

    #include "MpiCallbacks.hpp"
    void call_my_callback() {
      Communication::mpiCallbacks.call(my_callback, param1, param2);
    }

  This only works outside the integration loop. After the callback has been called, synchronous mpi communication can be done.

Legacy callbacks
~~~~~~~~~~~~~~~~

Older code uses callbacks defined in the CALLBACK_LIST preprocessor macro in communications.cpp. They are called via mpi_call().
See communications.cpp:mpi_place_particle() for an example.

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

The data structures for bonded interactions reside in ``interaction_data.hpp``.

* Add your interaction to the ``enum BondedInteraction``.
  This enumeration is used to identify different bonded interactions.
* Add a typedef struct containing the parameters of the interaction. Use the one for the FENE interaction as template::

    typedef struct {
      double k;
      [...]
    } Fene_bond_parameters;

* Add a member to the typedef union Bond_parameters. For the FENE bond it looks like this::

    Fene_bond_parameters fene;


Functions for calculating force and energy, and for setting parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every interaction resides in its own source .cpp and .hpp. A simple example for a
bonded interaction is the FENE bond in ``src/core/fene.cpp``` and ``src/core/fene.hpp``. 
Use these two files as templates for your interaction.

Notes:

* The names of function arguments mentioned below are taken from the FENE bond in ``src/core/feine.cpp`` and ``src/core/fene.hpp``. It is recommended to use the same names for the corresponding functions for your interaction. 
* The recommended signatures of the force and energy functions are::

    inline int calc_fene_pair_force(Particle *p1, Particle *p2, 
                                Bonded_ia_parameters *iaparams, 
                                double dx[3], double force[3])
    inline int fene_pair_energy(Particle *p1, Particle *p2, 
                            Bonded_ia_parameters *iaparams, 
                            double dx[3], double *_energy)

  Here, ``fene`` needs to be replaced by the name of the new interaction.
* The setter function gets a ``bond_type`` which is a numerical id identifying the number of the bond type in the simulation. It DOES NOT determine the type of the bond potential (harmonic vs FENE).
  The signature of the setter function has to contain the ``bond_type``, the remaining parameters are specific to the interaction. For the FENE bond, e.g., we have::

    fene_set_params(int bond_type, double k, double drmax, double r0)

  A return value of ``ES_OK`` is returned on success, ``ES_ERR`` on error, e.g., when parameters are invalid.
* The setter function must call make_bond_type_exists() with that bond type, to allocate the memory for storing the parameters.
* Afterwards, the bond parameters can be stored in the global variable bonded_ia_params[bond_type]
  
  * bonded_ia_params[bond_type].num is the number of particles involved in the bond -1. I.e., 1 for a pairwise bonded potential such as the FENE bond.
  * The parameters for the individual bonded interaction go to the member of Bond_parameters for your interaction defined in the previous step. For the FENE bond, this would be::
    bonded_ia_params[bond_tpe].p.fene
* At the end of the parameter setter function, do not forget the call to mpi_bcast_ia_params(), which will sync the parameters just set to other compute nodes in a parallel simulation.
* The routines for calculating force and energy return an integer. A return value of 0 means OK, a value of 1 means that the particles are too far apart and the bond is broken. This will stop the integration with a runtime error.
* The functions for calculating force and energy can make use of a pre-calculated distance vector (dx) pointing from particle 2 to particle 1.
* The force on particle 1 has to be stored in the force vector  (not added to it). The force on particle 2 will be obtained from Newton's law.
* The result of the energy calculation is placed in (NOT added to) the ``_energy`` argument of the energy calculation function.



Including the bonded interaction in the force calculation and the energy and pressure analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* In ``src/core/interaction_data.cpp``:

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

* Force calculation: in ``forces_inline.hpp`` in the function
  ``add_bonded_force()``, add your bond to the switch statement. For the FENE
  bond, e.g., the code looks like this::

    case BONDED_IA_FENE:
      bond_broken = calc_fene_pair_force(p1, p2, iaparams, dx, force);

* Energy calculation: add similar code to ``add_bonded_energy()`` in ``energy_inline.hpp``
* Pressure, stress tensor and virial calculation: If your bonded interaction is
  a pair bond and does not modify the particles involved, add similar code as
  above to pressure.hpp:calc_bonded_pair_force(). Otherwise, you have to
  implement a custom solution for virial calculation.


Adding the bonded interaciton in the Python interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please note that the following is Cython code (www.cython.org), rather than pure Python.
* In ``src/python/espressomd/interactions.pxd``:

  * import the parameter data structure from the C++ header file for your interaction. For the FENE bond, this looks like::

      cdef extern from "interaction_data.hpp":
          ctypedef struct Fene_bond_parameters:
              double k
              double drmax
              double r0
              double drmax2
              double drmax2i

  * Add your bonded interaction to the Cython copy of the BondedInteractions enum analogous to the one in the core:, described above::

      cdef enum enum_bonded_interaction "BondedInteraction":
          BONDED_IA_NONE = -1,
          BONDED_IA_FENE,
          BONDED_IA_HARMONIC,
          [...]

    The spelling has to match the one in the c++ enum exactly.
  * Adapt the Cython copy of the bond_parameters union analogous to the C++ core.  The member name has to match the one in C++ exactly::
      ctypedef union bond_parameters "Bond_parameters":
          Fene_bond_parameters fene
          Oif_global_forces_bond_parameters oif_global_forces
          Oif_local_forces_bond_parameters oif_local_forces
          Harmonic_bond_parameters harmonic
  * Import the declaration of the setter function implemented in the core. For the FENE bond, this looks like::
        cdef extern from "fene.hpp":
            int fene_set_params(int bond_type, double k, double drmax, double r0)

* In ``src/python/espressomd/interactions.pyx``:

  * Implement the Cython class for the bonded interaction, using the one for
    the FENE bond as template. Please use pep8 naming convention::

        class FeneBond(BondedInteraction):
        
            def __init__(self, *args, **kwargs):
                """ 
                FeneBond initialiser. Used to instatiate a FeneBond identifier
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
    
* In ``testsuite/python/bondedInteractions.py``:
  
  * Add a test case, which verifies that parameters set and gotten from the interaction are consistent::

        test_fene = generateTestForBondParams(
            0, FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})

  
  
  
   



.. _Outdated: Adding New Nonbonded Interactions:

Outdated: Adding New Nonbonded Interactions 
-------------------------------------------

Writing nonbonded interactions is similar to writing nonbonded
interactions. Again we start with ``interaction_data.h``, where the
parameter structure has to be set up. Just add your parameters *with
reasonable names* to ``IA_parameters``. Note that there must be a
setting for the parameters which disables the interaction.

Now write the header file for the interaction. This time ``ljcos.h`` may
be a good example. The needed routines are

-  ::

       int print*IAToResult(Tcl_Interp *interp, int i, int j)

   writes out the interaction parameters between particles of type ``i``
   and ``j`` to the interpreters result such that the result can be fed
   into the ``inter`` command again to obtain the same interaction. The
   ``IA_parameters`` pointer can be obtained conveniently via
   ``get_ia_param(i,j)``.

-  ::

       int *_parser(Tcl_Interp * interp, int part_type_a, int part_type_b, 
                    int argc, char ** argv)

   parses the command line given by ``argc`` and ``argv`` for the
   parameters needed for the interaction, and writes them to the
   ``IA_parameters`` for types ``part_type_a`` and ``part_type_b``. For
   details on writing the parser, see below. The routine returns 0 on
   errors and otherwise the number of parameters that were read from the
   command line.

-  ::

       void add_*_pair_force(Particle *p1, Particle *p2, 
                             IA_parameters *ia_params, 
                             double d[3], double dist2, double dist, 
                             double force[3])
       double *_pair_energy(Particle *p1, Particle *p2, 
                            IA_parameters *ia_params, 
                            double d[3], double dist2, double dist)

   are the routines to compute the force respectively the energy.
   ``ia_params`` gives the interaction parameters for the particle types
   of particles ``p1`` and ``p2``, ``d`` gives the vector from particle
   2 to particle 1, ``dist`` its length and ``dist2`` its squared
   length. The last three parameters can be chosen on demand. Note that
   unlike in the bonded case, the force routine is called ``add_*``,
   *i.e.*\ the force has to be *added* to force. The ``*_pair_energy``
   routine simply returns the energy directly instead of the pointer
   approach of the bonded interactions.

Change ``interaction_data.c`` as follows (most changes are pretty much
the same for all potentials):

#. modify ``initialize_ia_params`` and ``copy_ia_params`` to take care
   of the additional parameters needed for your potential.

#. ``checkIfParticlesInteract`` has to be modified to also check for the
   no interaction condition for the new interaction (typically zero
   cutoff).

#. ``calc_maximal_cutoff`` has to modified such that ``max_cut`` is
   larger than the maximal cutoff your interaction needs. Again, the
   code always uses the result from this function, therefore the cutoff
   does not have to be stored explicitly in the interaction parameters.

#. add your ``print*IAToResult`` routine to
   ``tclprint_to_result_NonbondedIA``.

#. add the ``*_parser`` routine to ``tclcommand_inter_parse_bonded``.

After this, add the force calculation to ``add_non_bonded_pair_force``,
``add_non_bonded_pair_virials`` and ``pressure_calc``, and the energy
calculation to ``add_non_bonded_pair_energy``.

After the new non-bonded interaction works properly, it would be a good
idea to add a testcase to the testsuite, so that changes breaking your
interaction can be detected early.

Outdated: Particle Data Organization
------------------------------------

The particle data organization is described in the Tcl command
cellsystem, its implementation is briefly described in ``cells.h`` and
``ghosts.h``. Here only some details on how to access the data is
assembled. Writing a new cellsystem almost always requires deep
interactions with the most low level parts of the code and cannot be
explained in detail here.

Typically, one has to access all real particles stored on this node, or
all ghosts. This is done via a loop similar to the following:

::

       Cell *cell;
       int c,i,np,cnt=0;
       Particle *part;
     
       for (c = 0; c < local_cells.n; c++) {
         cell = local_cells.cell[c];
         part = cell->part;
         np   = cell->n;
         for(i=0 ; i < np; i++) {
            do_something_with_particle(part[i]);
         }
       }

To access the ghosts instead of the real particles, use ``ghost_cells``
instead of ``local_cells``.

Another way to access particle data is via ``local_particles``. This
array has as index the particle identity, so that
``local_particles[25]`` will give you an pointer to the particle with
identity 25, or ``NULL``, if the particle is not stored on this node,
neither as ghost nor as real particle. Note that the ``local_particle``
array does not discriminate between ghosts and real particles. Its
primary use is for the calculation of the bonded interactions, where it
is used to efficiently determine the addresses of the bonding
partner(s).

The master node can add and remove particles via ``place_particle`` and
``remove_particle``, or change properties via ``set_particle_v`` etc.
This is the preferred way to handle particles, since it is
multiprocessor save.

However, some algorithms, especially new cellsystems, may force you to
operate locally on the particle data and shift them around manually.
Since the particle organization is pretty complex, there are additional
routines to move around particles between particle lists. The routines
exist in two versions, one indexed, and one unindexed. The indexed
version take care of the ``local_particles`` array, which for each
particle index tells where to find the particle on this node (or
``NULL`` if the particle is not stored on this node), while the
unindexed versions require you to take care of that yourself (for
example by calling ``update_local_particles``). The second way is much
faster if you do a lot of particle shifting. To move particles locally
from one cell to another, use ``move_indexed_particle`` or
``move_unindexed_particle``, never try to change something directly in
the lists, you will create a mess! Inserting particles locally is done
via ``append_indexed_particle`` or ``append_unindexed_particle``.

Besides the ``local_particles array``, which has to be up to date at any
time, there is a second array ``particle_node``, which is available on
the master node only outside of the integrator, *i.e.*\ in the Tcl
script evaluation phases. If ``particle_node`` is ``NULL``, you have to
call ``build_particle_node`` to rebuild it. For each particle identity
it contains the node that the particle is currently located on.

The proper cell for a particle is obtained via
``CellStructure::position_to_node``, which calculates for a given
position the node it belongs to, and
``CellStructure::position_to_cell``, which calculates the cell it
belongs to on this node, or ``NULL``, if the cell is from a different
node. However, you should normally not be bothered with this
information, as long as you stick to ``place_particle`` and the other
routines to modify particle data.

Writing a new cellsystem basically requires only to create the functions
listed in ``CellStructure``. The ``init`` function has to also setup the
communicators, which is the most complex part of writing a new
cellsystem and contains all the communication details. ``prepare_comm``
is a small wrapper for the most common operations. Otherwise just grep
for ``CELL_STRUCTURE_DOMDEC``, and add some appropriate code for your
cell system. Note, however, that each cell system has its specific part
of the code, where only this cellsystem does something strange and
unique, so here you are completely on your own. Good luck.

.. _Outdated\: Errorhandling for Developers:

Outdated: Errorhandling for Developers
--------------------------------------

Developers should use the errorhandling mechanism whenever it is
possible to recover from an error such that continuing the simulation is
possible once the source of the error is removed, i. e. the bond is
removed or a parameter changed. For example, if due to excessive forces,
particles have been far out of their current node, |es| puts them into
one of the local cells. Since the position is unphysical anyways, it is
of no importance anymore, but now the user can place the particles anew
and perhaps decrease the time step such that the simulation can continue
without error. However, most often the recovery requires no special
action.

To issue a background error, call

::

    errtxt=runtime_error(length)

where length should be the maximal length of the error message (you can
use ``TCL_DOUBLE_SPACE`` rsp. ``TCL_INTEGER_SPACE`` to obtain space for
a double rsp. integer). The function returns a pointer to the current
end of the string in ``error_msg``. After doing so, you should use the
``ERROR_SPRINTF``-macro, which substitutes to a simple ``sprintf``, so
that your errormessage will automatically be added to the
“runtime-errors resolved”-page. Please make sure that you give each of
your errors an unique 3-digit errorcode (for already used errorcodes
have a look at the “runtime-errors resolved”-page), have the curled
braces around your message and the space at the end, otherwise the final
error message will look awful and will propably not automatically be
added to our error-page. Typically, this looks like this::

    if (some_error_code != OK) {
      char *errtxt = runtime_error(TCL_INTEGER_SPACE + 128);
      ERROR_SPRINTF(errtxt, "{error occured %d} ", some_error_code);
      recovery;
    }

If you have long loops during which runtime errors can occur, such as
the integrator loop, you should call ``check_runtime_errors`` from time
to time and exit the loop on errors. Note that this function requires
all nodes to call it synchronously.

In all cases, all Tcl commands should call ``mpi_gather_runtime_errors``
before exiting. You simply handover the result you were just about to
return. If the result was ``TCL_ERROR``, then
``mpi_gather_runtime_errors`` will keep the Tcl error message and
eventually append the background errors. If the result was ``TCL_OK``,
*i.e.*\ your function did not find an error, the result will be reset
(since |es| is in an undefined state, the result is meaningless), and
only the background errors are returned. Whenever a Tcl command returns,
instead of ``return TCL_OK/TCL_ERROR`` you should use

::

    return mpi_gather_runtime_errors(interp, TCL_OK/TCL_ERROR); 






Global Variables which are synchronized across nodes
----------------------------------------------------

Adding new global variables to |es|, is strongly discouraged, because it means that code depends on a purely defined global state and cannot be tested individually.
Features/Algorithms should instead be encapsulated in a class which is used by the script interface mechanism.

However, there is a mechanism in the simulation core, to synchronize existing global variables across the mpi cores.

These variables are declared ``extern`` in a header file and include in
``global.cpp``. Then there is a line to the definition of the constant data
structure ``fields`` at the beginning of the file ``global.c``. For
details on the entries, see the definition of ``Datafield`` in
``global.h``). Basically it is declare *where* the variable is
stored, *which type* (INT or DOUBLE) it has and *how many* elements. A
callback procedure can be provided which checks if the given value is
valid and stores it. It is also responsible for dispatching the new
value to the other compute nodes, if necessary. This is done via ``mpi_bcast_parameter()``, which will transfer the value
to the other nodes. A simple example is ``box_l`` with the callback
procedure ``boxl_callback``. For ``mpi_bcast_parameter`` to work, it is
necessary that they occur in the list of constant definitions at the
beginning of ``global.hpp``. So please keep this list in sync!


.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
