=================
Developer's Guide
=================

.. warning::
   The information found in this version of the Developer's Guide is
   outdated.  Please see the section :ref:contact` and
   ask for advice if you plan to start working on |es|.

.. _contact:

Contact the Develoeprs
======================

To contact the |es| developers, please write an email to the developerss mailing list:
espressomd-devel@nongnu.org
to subscribe to the developers’ mailing list go to
http://lists.nongnu.org/mailman/listinfo/espressomd-devel 


Before you start a development project
======================================
Before you start a development project for |es|, please always write to the developers mailing list and describe the project. 
This is to avoid that several people work on the same thing at the same time. Also, implementation details can be discussed in advance. In many cases, existing developers can point to re-usable code and simpler solutions.



.. _development_environment:

Development Environment
=======================

.. _required_development_tools:

Required Development Tools
--------------------------
-  First of all, please install the dependencies for compiling |es|. See :ref:`_Getting, comping and running`

-  To be able to access the development version of |es|, you will need
   the distributed versioning control system Git [1]_. Section
   :ref:`git_repositories` contains documentation on how we employ
   git.

-  To be able to compile the User’s Guide or the Developer’s Guide, you
   will need a LaTeX-installation (all recent ones will do). For
   details, refer to chapter .

-  To compile the Doxygen code documentation, you will need to have the
   tool Doxygen\  [4]_.

All of these tools should be easy to install on most Unix operating
systems.

.. _getting_the_development_code:

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

Adding New SOURCE Files
^^^^^^^^^^^^^^^^^^^^^^^

To add new files to |es| (like C++ source files or header files) you
need to look at the CMakeList.txt in the directory where the file is located.
-  In some cases (e.g., src/core/CMakeList.txt), the CMakeList.txt contains a wild-card include like this
   file(GLOB EspressoCore_SRC
          "*.cpp"
          )
   In this case, placing a file with that ending is enough.
-  Please note that .hpp-header files usually do not have to be added to CMakeList.txt
-  In other cases, the files are explicitly included (e.g., testsuite/python/CMakeList), 
   set(py_tests  bondedInteractions.py
              cellsystem.py
              constraint_shape_based.py
              coulomb_cloud_wall.py
   In that case, add the new file to the list.
   


Testsuite
---------
-  New or significantly changed features will only be accepted, if they have a test case. 
   This is to make sure, the feature is not broken by future changes to |es|, and so other users can get an impression of what behaviour is guaranteed to work.
-  There are two kinds of tests:

  -  C++-unit tests, testing individual c++ functions and classes. They make use of the boost unit test framework and reside in ``src/core/unit_tests`
  -  Python integration tests, testing the Python interface and (physical) results of features. They reside in ``testsuite/python``

-  To execute the tests, run
   make check 
   in the top build directory.


.. _documentation:

Documentation
=============
The documentation of |es| consists of three parts:
-  The users' guide and developers' guide are located in ``doc/sphinx``, and make use of the Sphinx Python package
-  In-code documentation for the Python interface is located in the various files in src/python/espressomd and also makes use of the Sphinx Python package
-  In-code documentation of the C++ core is located in the .cpp and .hpp files in ``/sr/core`` and its sub-directories and makes use of Doxygen.



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
     pictures, as doxygen will not copy the pictures into the
     documentation.

-  | ``<ul> <li>List entry 1</li> <li>List entry 2</li></ul>``
   | Creates a list in the documentation.

-  | ``\param`` *name* *description*
   | Document the parameter of a function.

-  | ``\return`` *decription*
   | Document the return value of a function.

.. _programmers_guide:

Programmer’s Guide
==================

This chapter provides some hints on how to extend |es|. It is not
exhaustive, so for major changes the best documentation are the other
developers.

.. _adding_global_variables:

.. _adding_new_bonded_interactions:

Adding New Bonded Interactions
------------------------------

Every interaction resides in its own source file. A simple example for a
bonded interaction is the FENE bond in ``fene.h``. The data structures,
however, reside in ``interaction_data.h``. The bonded interactions are
all stored in a union, ``Bonded_ia_parameters``. For a new interaction,
just add another struct. Each bonded interaction is assigned a type
number, which has the form ``BONDED_IA_*``, *e.g.*\ ``BONDED_IA_FENE``.
The new interaction also has to have such a *unique* number.

After the setup of the necessary data structures in
``interaction_data.h``, write the source file, something like
``new_interaction.h``. You may want to use ``fene.h`` as a template
file. Typically, you will have to define the following procedures:

-  ::

       int *_set_params(int bond_type, ...)

   This function is used to define the parameters of a bonded
   interaction. ``bond_type`` is the bond type number from the inter
   command, and not one of the ``BONDED_*``. It is rather an index to
   the ``bonded_ia_params`` array. ``make_bond_type_exist`` makes sure
   that all bond types up the given type exist and are preinitialized
   with ``BONDED_IA_NONE``, *i.e.*\ are empty bond types. Therefore fill
   ``bonded_ia_params[bond_type]`` with the parameters for your
   interaction type.

-  ::

       int calc_*_force(Particle *p1, Particle *p2,..., 
                        Bonded_ia_parameters *iaparams, 
                        double dx[3], double force[3], ...)

   This routine calculate the force between the particles. ``ia_params``
   represents the parameters to use for this bond, ``dx`` represents the
   vector pointing from particle 2 to particle 1. The force on particle
   1 is placed in the force vector (and *not* added to it). The force on
   particle 2 is obtained from Newton’s law. For many body interactions,
   just add more particles in the beginning, and return the forces on
   particles 1 to N-1. Again the force on particle N is obtained from
   Newton’s law. The procedure should return 0 except when the bond is
   broken, in which case 1 is returned.

-  ::

       int *_energy(Particle *p1, Particle *p2, ..., 
                    Bonded_ia_parameters *iaparams, 
                    double dx[3], double *_energy)

   This calculates the energy originating from this bond. The result is
   placed in the location ``_energy`` points to, ``ia_params`` and
   ``dx`` are the same as for the force calculation, and the return
   value is also the flag for a broken bond.

After the preparation of the header file, the bonded interaction has to
be linked with the rest of the code. In ``interaction_data.c``, most of
the work has to be done:

#. Add a name for the interaction to ``get_name_of_bonded_ia``.

#. In ``calc_maximal_cutoff``, add a case for the new interaction which
   makes sure that ``max_cut`` is larger than the interaction range of
   the new interaction, typically the bond length. This value is always
   used as calculated by ``calc_maximal_cutoff``, therefore it is not
   strictly necessary that the maximal interaction range is stored
   explicitly.

#. Add a print block for the new interaction to
   ``tclcommand_inter_print_bonded``. The print format should be such
   that the output can be used as input to inter, and defines the same
   bond type.

#. In ``tclcommand_inter_parse_bonded``, add a parser for the
   parameters. See the section on parsing below.

#. Besides this, you have enter the force respectively the energy
   calculation routines in ``add_bonded_force``, ``add_bonded_energy``,
   ``add_bonded_virials`` and ``pressure_calc``. The pressure occurs
   twice, once for the parallelized isotropic pressure and once for the
   tensorial pressure calculation. For pair forces, the pressure is
   calculated using the virials, for many body interactions currently no
   pressure is calculated.

After the new bonded interaction works properly, it would be a good idea
to add a testcase to the testsuite, so that changes breaking your
interaction can be detected early.

.. _adding_new_nonbonded_interactions:

Adding New Nonbonded Interactions
---------------------------------

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

Particle Data Organization
--------------------------

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

.. _errorhandling_for_developers:

Errorhandling for Developers
----------------------------

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
added to our error-page. Typically, this looks like this:

::

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
-----------------------------------------------------------

Adding new global variables to |es|, is strongly discuraged, because it means that code depends on a purely defined global state and cannot be tested individually.
Features/Algorithms should instead be encapsulated in a class which is used by the script interface mechnaism.

However, there is a mechanism in the simulation core, to synchronize existing global variables across the mpi cores.

these variables are declared ``extern`` in a header file and include in
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


.. [1]
   http://git-scm.com/

.. [2]
   http://www.gnu.org/software/automake/

.. [3]
   http://www.gnu.org/software/autoconf/autoconf.html

.. [4]
   http://www.doxygen.org/

.. [5]
   http://www.gnu.org/software/automake/

.. [6]
   http://www.gnu.org/software/autoconf/autoconf.html
